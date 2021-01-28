function solver_logging(iter,α,S,dS)
    if rem((iter-1),4)==0
    printstyled("iter     α              S              dS              \n";
                    bold = true, color = :light_blue)
        bars = "----------------------------"
        bars*="---------------------\n"
        printstyled(bars; color = :light_blue)
    end
    # print(Crayon(foreground = :red), "In red. ", Crayon(bold = true), "Red and bold")
    DJ = @sprintf "%.4E" dS
    maxL = @sprintf "%.4E" S
    J_display = @sprintf "%.4E" dS
    alpha_display = @sprintf "%.4E" α
    str = "$iter"
    for i = 1:(6 - ndigits(iter))
        str *= " "
    end
    println("$str   $alpha_display     $maxL     $J_display     ")
    return nothing
end
# function solver_logging2(iter,DJ,l,J,α)
#     """Don't worry about this, just for visual logging in the REPL"""
#     if rem((iter-1),4)==0
#     printstyled("iter     α              maxL           J              DJ\n";
#                     bold = true, color = :light_blue)
#         bars = "----------------------------"
#         bars*="------------------------------------\n"
#         printstyled(bars; color = :light_blue)
#     end
#     DJ = @sprintf "%.4E" DJ
#     maxL = @sprintf "%.4E" round(maximum(maximum.(l)),sigdigits = 3)
#     J_display = @sprintf "%.4E" J
#     alpha_display = @sprintf "%.4E" α
#     str = "$iter"
#     for i = 1:(6 - ndigits(iter))
#         str *= " "
#     end
#     println("$str   $alpha_display     $maxL     $J_display     $DJ")
#     return nothing
# end
