function residual(eci_initial)
    eci_i = copy(eci_initial)
    eci_i[1:3] *= dscale
    eci_i[4:6] *= (dscale/tscale)
    eci_hist_cubic = rollout(eci_i)
    outp = zeros(M*3)
    for i = 1:M
        t = Y_ts[i]
        # @infiltrate
        # error()
        outp[idx_r[i]] = measurement(eci_hist_cubic(t),t) - Y[i]
    end
    return outp
end

function GN(v0)
    println(" ")
    println("---------------BEGIN GAUSS-NEWTON----------------")
    new_S = 1e15
    for i = 1:10

        res = residual(v0)
        S = dot(res,res)
        # @show S

        # print initial cost
        if i == 1
            S_p = @sprintf "%.6E" S
            printstyled("Initial cost: $S_p\n",bold = true, color = :magenta)
        end

        # jacobian
        J = FiniteDiff.finite_difference_jacobian(residual,v0)
        # newton_step = -J\res
        newton_step = -(J'*J + 1e-6*I)\(J'*res)
        α = 1.0
        dS = NaN
        # backtracking line search
        for ii = 1:20
            new_v = v0 + α*newton_step
            new_res = residual(new_v)
            new_S = dot(new_res,new_res)
            if new_S < S
                dS = abs(S - new_S)
                v0 = copy(new_v)
                break
            else
                α /=2
            end
            if ii == 20
                dS = 0
            end
        end
        res = residual(v0)
        S = dot(res,res)
        solver_logging(i,α,S,dS)
        if dS<1e-5
            break
        end
    end
    return v0
end
