################################
# Parameterize OTN Model
################################

using DataFrames, CSV, LinearAlgebra, Statistics, StatsBase, Plots, NaNStatistics
using OptimalTransportNetworks

# Read Undirected Graph
graph = CSV.read("data/trans_africa_network/trans_african/fastest_routes_graph_edges.csv", DataFrame)
# histogram(graph.ug_cost, bins=100)
# Adjusting 
graph.distance /= 1000    # Convert to km
graph.border_dist /= 1000 # Convert to km
# graph.sp_distance /= 1000 # Convert to km
graph.duration /= 60      # Convert to hours
graph.total_cost = graph.ug_cost / 1e6  # Convert to millions

n = maximum([maximum(graph.from), maximum(graph.to)])

# Create Adjacency Matrix
adj_matrix = falses(n, n)
for i in 1:size(graph, 1)
    adj_matrix[graph.from[i], graph.to[i]] = adj_matrix[graph.to[i], graph.from[i]] = true
end

# Read Nodes Data
nodes = CSV.read("data/trans_africa_network/trans_african/fastest_routes_graph_nodes.csv", DataFrame)
nodes.population /= 1000 # Convert to thousands
nodes.outflows /= 1000 # Convert to thousands
describe(nodes)

# Create Infrastructure Matrix: Following Graff (2024) = average speed in km/h: length of route is accounted for in cost function
infra_matrix = zeros(n, n)
for i in 1:size(graph, 1)
    speed = graph.distance[i] / graph.duration[i]
    infra_matrix[graph.from[i], graph.to[i]] = infra_matrix[graph.to[i], graph.from[i]] = speed
end
describe(graph.distance ./ graph.duration)

# Create Iceberg Trade Cost Matrix: Following Graff (2024)
# delta_0_tau_withcomp_fp <- (0.0248/(0.43*0.25541) + 0.0254/(1.03*4.1418)) / 2 # this is the mean of columns (2) and (5) on page 44 of their paper, scaled by mean base prices for each country, so as to make these ad-valorem
# delta_tau_withcomp_fp <- delta_0_tau_withcomp_fp * log(dist / 1.609)
# delta_tau_withcomp_fp[delta_tau_withcomp_fp < 0] <- 0
iceberg_matrix = zeros(n, n)
graph.distance .+= graph.border_dist
for i in 1:size(graph, 1)
    iceberg_matrix[graph.from[i], graph.to[i]] = iceberg_matrix[graph.to[i], graph.from[i]] = 0.1158826 * log(graph.distance[i] / 1.609)
end
iceberg_matrix[iceberg_matrix .< 0] .= 0

# From Collier et al. (2016):
infra_building_matrix = zeros(n, n)
for i in 1:size(graph, 1)
    infra_building_matrix[graph.from[i], graph.to[i]] = infra_building_matrix[graph.to[i], graph.from[i]] = graph.total_cost[i]
end

# Basic characteristics of the economy
population = nodes.population;
population += (population .== 0) * 1e-6;

# Productivity Matrix 
# Check largest pcities
sum(population .> 2000 .|| nodes.outflows .> 1000) == 47
productivity = zeros(n, maximum(nodes.product))
with_ports = true
for i in 1:n
    productivity[i, nodes.product[i]] = nodes.IWI[i]
    if with_ports && nodes.outflows[i] > 0
        productivity[i, nodes.product[i]] += (37 * nodes.outflows[i]) / population[i] # Productivity of ports
    end
end
extrema(productivity)
all(sum(productivity .> 0, dims = 2) .== 1)

J = size(productivity, 1);
N = size(productivity, 2);

## Optimisation
min_mask = infra_matrix;                          
max_mask = max.(infra_matrix, adj_matrix .* 100);   

K_inv = sum(infra_building_matrix) / 2 # Cost of all Road Activities (33 billion) 
# Total implied costs
infra_building_matrix ./= (max_mask - infra_matrix)
infra_building_matrix[isinf.(infra_building_matrix)] .= 0
infra_building_matrix[isnan.(infra_building_matrix)] .= 0
K_base = sum(infra_building_matrix .* infra_matrix) / 2
K = (K_base + 10e3) * 2 # 10 or 20 billion investment volume (*2 because symmetric)

K / sum(infra_building_matrix .* max_mask)
K / sum(infra_building_matrix .* min_mask)
if K > sum(infra_building_matrix .* max_mask)
    error("Infrastructure budget is too large for the network")
end

# Parameters
alpha = 0.7           # 0.4? # Curvature parameter of the utility function 
a = 1                 # TG: 0.7
# rho = 2 # With inequality aversion !! 
gamma_DRS = 0.946         # TG
beta = 1.2446 * gamma_DRS # TG
gamma_IRS = beta^2/gamma_DRS  # IRS

# TODO: Debug Annealing procedure: why does solver return invalid model??
for gamma in [gamma_IRS], rho in [0, 2], sigma in [2, 5]
    print("\n\n\n\nsigma = ", sigma, " gamma = ", gamma, " rho = ", rho, "\n\n")

    ## Initialise geography
    param = init_parameters(annealing = false, labor_mobility = false, cross_good_congestion = true, 
                            a = a, sigma = sigma, N = N, alpha = alpha, beta = beta, 
                            gamma = gamma, rho = rho, K = K, 
                            tol = 1e-5, min_iter = 20, max_iter = 60);

   g = create_graph(param, type = "custom", x = nodes.lon, y = nodes.lat, adjacency = adj_matrix, 
                    Lj = population, Zjn = productivity, Hj = population .* (1-alpha)); # I normalise this because the general utility function has a (h_j/(1-alpha))^(1-alpha) thing with it

    g[:delta_i] = infra_building_matrix;
    g[:delta_tau] = iceberg_matrix;

    param[:optimizer_attr] = Dict(:hsllib => "/usr/local/lib/libhsl.dylib", :linear_solver => "ma86")

    # Static
    @time res_stat = optimal_network(param, g, I0 = infra_matrix, verbose = true, solve_allocation = true); # I0 = infra_matrix, 

    # Optimal 
    @time res_opt = optimal_network(param, g, I0 = infra_matrix, Il = min_mask, Iu = max_mask, verbose = false); # I0 = infra_matrix, 

    # Check
    print(K / sum(res_opt[:Ijk] .* g[:delta_i]))

    # Plot Network
    display(plot_graph(g, res_opt[:Ijk], height = 800)) #  # , node_sizes = res[:Cj])
    display(plot_graph(g, res_opt[:Ijk] - infra_matrix, height = 800))

    # Save
    function res_to_vec(Ijk, graph)
        n = size(graph, 1)
        rv = zeros(n)
        for i in 1:n
            rv[i] = (Ijk[graph.from[i], graph.to[i]] + Ijk[graph.to[i], graph.from[i]]) / 2
        end
        return rv
    end

    # Saving: Nodes
    res_nodes = deepcopy(nodes)
    res_nodes.uj_orig = res_stat[:uj]
    res_nodes.Lj_orig = res_stat[:Lj]
    res_nodes.Cj_orig = res_stat[:Cj]
    res_nodes.Dj_orig = res_stat[:Dj]
    res_nodes.PCj_orig = res_stat[:PCj]
    # res_nodes.welfare = res_opt.welfare; sum(res_opt.Lj .* res_opt.uj)
    res_nodes.uj = res_opt[:uj]
    res_nodes.Lj = res_opt[:Lj]
    res_nodes.Cj = res_opt[:Cj]
    res_nodes.Dj = res_opt[:Dj]
    res_nodes.PCj = res_opt[:PCj]
    for n in 1:N
        res_nodes[!, Symbol("Lj_$(n)")] = res_opt[:Ljn][:,n]
        res_nodes[!, Symbol("Dj_$(n)")] = res_opt[:Djn][:,n]
        # res_nodes[!, Symbol("Cj_$(n)")] = res_opt[:Cjn][:,n]
        res_nodes[!, Symbol("Yj_$(n)")] = res_opt[:Yjn][:,n]
        res_nodes[!, Symbol("Pj_$(n)")] = res_opt[:Pjn][:,n]
    end
    res_nodes |> CSV.write("results/trans_african/nodes_results_22g$(with_ports ? "" : "_noport")_10b_fixed_cgc$(gamma == gamma_IRS ? "_irs_na" : "")_sigma$(sigma)_rho$(rho)_bc_julia.csv")

    # Saving: Graph
    res_graph = deepcopy(graph)
    res_graph.Ijk_orig = res_to_vec(infra_matrix, graph)
    res_graph.Ijk = res_to_vec(res_opt[:Ijk], graph)
    for n in 1:N
        res_graph[!, Symbol("Qjk_$(n)")] = res_to_vec(res_opt[:Qjkn][:,:,n], graph)
    end
    res_graph |> CSV.write("results/trans_african/edges_results_22g$(with_ports ? "" : "_noport")_10b_fixed_cgc$(gamma == gamma_IRS ? "_irs_na" : "")_sigma$(sigma)_rho$(rho)_bc_julia.csv")
end



