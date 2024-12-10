using DelimitedFiles

mutable struct Parameters
    name::String
    alpha::Float64
    beta::Float64
    gamma::Float64
    eta::Float64
    lambda::Float64
    mu::Float64
    T_0::Float64
    T_S::Float64
    K_N::Float64
    K_S::Float64
    K_IP::Float64
    S_0::Float64
    S_N::Float64
    S_T::Float64
    S_S::Float64
    S_IP::Float64
    S_B::Float64
    F_N::Float64
    F_T::Float64
    F_S::Float64
    F_IP::Float64
    A_N::Float64
    A_T::Float64
    A_S::Float64
    A_IP::Float64
    Y::Float64
    C::Float64
    V_N::Float64
    V_T::Float64
    V_S::Float64
    V_IP::Float64
    V_B::Float64
    on_state::Vector{Float64}
    off_state::Vector{Float64}
    q_on::Float64
    q_off::Float64
    on_state_5::Vector{Float64}
    off_state_5::Vector{Float64}
    q_on_5::Float64
    q_off_5::Float64
    B11::Float64
    B21::Float64
    B22::Float64
    noise_type::String
    noise_scaling::Int64
end

function create_param(paramset)
    current_dir = @__DIR__
    filename = joinpath(current_dir, "parameters.txt")
    parameters = readdlm(filename)

    # Extract parameters from file
    V_N = parameters[1, paramset] * 1e16
    V_T = parameters[2, paramset] * 1e16
    V_S = parameters[3, paramset] * 1e16
    V_IP = parameters[4, paramset] * 1e16
    V_B = parameters[5, paramset] * 1e16
    T_S = parameters[14, paramset]
    T_0 = parameters[15, paramset]
    lambda = parameters[17, paramset] * 1e7
    gamma = parameters[22, paramset]
    mu = parameters[16, paramset] * 1e-8
    eta = parameters[21, paramset] * 1e6
    K_N = parameters[18, paramset] * 1e6
    K_S = parameters[19, paramset] * 1e6
    K_IP = parameters[20, paramset] * 1e6
    F_N = parameters[10, paramset] * 1e6
    F_T = parameters[12, paramset] * 1e6
    F_S = parameters[11, paramset] * 1e6
    F_IP = parameters[13, paramset] * 1e6
    A_N = parameters[6, paramset] * 1e6
    A_T = parameters[7, paramset] * 1e6
    A_S = parameters[8, paramset] * 1e6
    A_IP = parameters[9, paramset] * 1e6


    # Define fix parameters
    alpha = 0.12 # kg m^-3 C^-1
    beta = 790.0 # kg m^-3 
    Y = 3.15e7 # sec/year
    S_N_eq = 0.034912
    S_T_eq = 0.035435
    S_S_eq = 0.034427
    S_IP_eq = 0.034668
    S_B_eq = 0.034538

    C = V_N * S_N_eq + V_T * S_T_eq + V_S * S_S_eq + V_IP * S_IP_eq + V_B * S_B_eq

    S_0 = 0.035
    S_N = 100*(S_N_eq-S_0)
    S_T = 100*(S_T_eq-S_0)
    S_S = 100*(S_S_eq-S_0)
    S_IP = 100*(S_IP_eq-S_0)
    S_B = 100*(S_B_eq-S_0)

    if paramset == 1 # FamousA 1×CO2
        name = "FAMOUS_A 1×CO2"
        on_state = [0.0071673, 0.065394]
        off_state = [NaN, NaN] # [-0.117257, 0.0980277] # no equillibrium but a valley
        q_on = 1.710248e7
        q_off = NaN
        on_state_5 = [NaN, NaN, NaN, NaN]
        off_state_5 = [NaN, NaN, NaN, NaN]
        q_on_5 = NaN
        q_off_5 = NaN
    elseif  paramset == 2 # FamousB 1×CO2
        name = "FAMOUS_B 1×CO2"
        on_state = [-0.0082838, 0.0570158]
        off_state = [-0.1180082, 0.0586036]
        q_on = 1.512614e7
        q_off = -5.29742e6
        on_state_5 = [-0.0056424, 0.0583561, -0.0569066, -0.0317845]
        off_state_5 = [-0.1151225, 0.0526256, -0.0488454, -0.0192607]
        q_on_5 = 1.554458e7
        q_off_5 = -6.33399e6
    elseif paramset == 3 # FamousB 2×CO2
        name = "FAMOUS_B 2×CO2"
        on_state = [0.0323752, 0.1434514]
        off_state = [-0.1984548, 0.1499503]
        q_on = 1.355190e7
        q_off = -7.14007e6
        on_state_5 = [0.027207, 0.1270648, -0.0816578, -0.0545624]
        off_state_5 = [-0.1958341, 0.1440135, -0.0514548, -0.017383]
        q_on_5 = 1.5272170e7
        q_off_5 = -7.42912e6
    elseif paramset == 4 # HadGEM2-AO 1×CO2
        name = "HadGEM2-AO 1×CO2"
        on_state = [-0.0253798, 0.0527541]
        off_state = [-0.1456736, 0.0449482]
        q_on = 1.469108e7
        q_off = -5.20554e6
        on_state_5 = [NaN, NaN, NaN, NaN]
        off_state_5 = [NaN, NaN, NaN, NaN]
        q_on_5 = NaN
        q_off_5 = NaN
    elseif paramset == 5 # HadGEM2-AO 2×CO2
        name = "HadGEM2-AO 2×CO2"
        on_state = [0.0013263, 0.0634038]
        off_state = [NaN, NaN]
        q_on = 1.207473e7
        q_off = NaN
        on_state_5 = [NaN, NaN, NaN, NaN]
        off_state_5 = [NaN, NaN, NaN, NaN]
        q_on_5 = NaN
        q_off_5 = NaN
    elseif paramset == 6 # HadGEM2-AO 4×CO2
        name = "HadGEM2-AO 4×CO2"
        on_state = [0.0129943, 0.0775040]
        off_state = [NaN, NaN]
        q_on = 1.035154e7
        q_off = NaN
        on_state_5 = [NaN, NaN, NaN, NaN]
        off_state_5 = [NaN, NaN, NaN, NaN]
        q_on_5 = NaN
        q_off_5 = NaN
    elseif paramset == 7 # HadGEM3LL
        name = "HadGEM3-LL"
        on_state = [-0.0005551, 0.1216170]
        off_state = [-0.1318978, 0.0333245]
        q_on = 8.592036e6
        q_off = -8.943528e6
        on_state_5 = [NaN, NaN, NaN, NaN]
        off_state_5 = [NaN, NaN, NaN, NaN]
        q_on_5 = NaN
        q_off_5 = NaN
    elseif paramset == 8 # HadGEM3MM
        name = "HadGEM3-MM"
        on_state = [-0.0339449, 0.1048730] # no stable on_state
        off_state = [-0.1496008, -0.0260129]
        q_on = 6.627932e6
        q_off = -1.464252e7
        on_state_5 = [NaN, NaN, NaN, NaN]
        off_state_5 = [NaN, NaN, NaN, NaN]
        q_on_5 = NaN
        q_off_5 = NaN
    else
        name = ""
        on_state = [NaN, NaN]
        off_state = [NaN, NaN]
        q_on = NaN
        q_off = NaN
        on_state_5 = [NaN, NaN, NaN, NaN]
        off_state_5 = [NaN, NaN, NaN, NaN]
        q_on_5 = NaN
        q_off_5 = NaN
    end

    noise_factor = 10

    # HadGEM3MM
    noise_type = "HadGEM3-MM"
    B11 =  0.1263 * 1e-3 * noise_factor
    B21 = -0.0869 * 1e-3 * noise_factor
    B22 =  0.1088 * 1e-3 * noise_factor

    # Create a new instance of Parameters
    return Parameters(name, alpha, beta, gamma, eta, lambda, mu, T_0, T_S, K_N, K_S, K_IP, S_0, S_N, S_T, S_S, S_IP, S_B, F_N, F_T, F_S, F_IP, A_N, A_T, A_S, A_IP, Y, C, V_N, V_T, V_S, V_IP, V_B, on_state, off_state, q_on, q_off, on_state_5, off_state_5, q_on_5, q_off_5, B11, B21, B22, noise_type, noise_factor)

end

mutable struct lorenz_parameters
    rho::Float64
    mu::Float64
    beta::Float64
    epsilon::Float64 # sqrt(vel^-1)
    delta::Float64 # strength = sqrt(2δ)
    sigma::Float64
end

mutable struct forcing_parameters
    forcing::Bool
    initial_forcing::Float64
    final_forcing::Float64
    height_forcing::Float64
    t_f::Float64
    t_rise::Float64
    t_pert::Float64
    t_fall::Float64
end
