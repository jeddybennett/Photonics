using QuantumOptics
using Plots 
using BenchmarkTools
using LinearAlgebra
##
# bases, states, tensor products
##

b2 = SpinBasis(1//2)
bf = FockBasis(5)

su = spinup(b2)
f1 = fockstate(bf, 1)

ψᵢ = su ⊗ f1

n = number(bf)
a = destroy(bf)
ad = a'
id_f = identityoperator(bf)

σx = sigmax(b2)
σz = sigmaz(b2)
id_2 = identityoperator(b2)

##
# partial traces, embedding
##

ρᵢ = dm(ψᵢ)

ptrace(ρᵢ, 1)

embed(basis(ψᵢ), [2], n) == id_2⊗n


function prepare_test(n)
    b = FockBasis(n)
    ϕ = coherentstate(b, 1.0)
    op = create(b)
    return ϕ, op
end

function test_performance(op, state)
    return op*state
end

n_for_tests = [2, 3, 10, 100, 500]
benchmark_results_sparse = []
for n in n_for_tests
    println(n)
    ϕ, op = prepare_test(n)
    println(op)
    b = @benchmark test_performance($op, $ϕ)
    t = mean(b).time
    push!(benchmark_results_sparse, t)
end

benchmark_results_dense = []
for n in n_for_tests
    println(n)
    ϕ, op = prepare_test(n)
    op = dense(op)
    b = @benchmark test_performance($op, $ϕ)
    t = mean(b).time
    push!(benchmark_results_dense, t)
end

#schordinger equation
##

tspan = 0:100.0

plusHC(op) = op + op'

H = id_2 ⊗ n + plusHC(sigmam(b2)⊗create(bf))

out_t, out_states = timeevolution.schroedinger(tspan, ψᵢ, H)

expectation_number = [abs(expect(2, n, s)) for s in out_states]

plot(out_t, expectation_number)

##

#lindblad master equation
##


tspan = 0:100.0

plusHC(op) = op + op'

H = id_2 ⊗ n + plusHC(sigmam(b2)⊗create(bf))

#randomize phase of your spin space
J = [0.1*sigmaz(b2)⊗id_f, 0.1*id_2⊗destroy(bf)]

out_t, out_states = timeevolution.master(tspan, ψᵢ, H, J)

expectation_number = [abs(expect(2, n, s)) for s in out_states]

#meaning our photon number is decreasing
plot(out_t, expectation_number)


##
#monte carlo trajectories
tspan = 0:100.0

plusHC(op) = op + op'

H = id_2 ⊗ n + plusHC(sigmam(b2)⊗create(bf))

#randomize phase of your spin space
J = [0.1*sigmaz(b2)⊗id_f, 0.1*id_2⊗destroy(bf)]
expectation_number_repeated_samples = []
for sample in 1:100
    out_t, out_states = timeevolution.mcwf(tspan, ψᵢ, H, J)
    expectation_number = [abs(expect(2, n, s)) for s in out_states]
    push!(expectation_number_repeated_samples, expectation_number)
end
#trajectories of sampling

##


plot()
for ex_n_samp in expectation_number_repeated_samples
    plot!(out_t, ex_n_samp, label=nothing, alpha=0.1)
end


##


avg_traj = sum(expectation_number_repeated_samples)/length(expectation_number_repeated_samples)
plot(out_t, avg_traj)


##


# time dependent operators 
tspan = 0:100.0

plusHC(op) = op + op'

function Htimedep(t, psi)
    id_2 ⊗ n + sin(t)*plusHC(sigmam(b2)⊗create(bf))
end

##

out_t, out_states = timeevolution.schroedinger_dynamic(tspan, ψᵢ, Htimedep)

expectation_number = [abs(expect(2, n, s)) for s in out_states]

plot(out_t, expectation_number)

##

@time out_t, out_states = timeevolution.schroedinger_dynamic(tspan, ψᵢ, Htimedep);


##

Hlazytimedep = TimeDependentSum(
    [1.0, sin],
    [id_2 ⊗ n, plusHC(sigmam(b2)⊗create(bf))]
)


@time out_t, out_states = timeevolution.schroedinger_dynamic(tspan, ψᵢ, Htimedep);
@time out_t, out_states = timeevolution.schroedinger_dynamic(tspan, ψᵢ, Hlazytimedep);



##