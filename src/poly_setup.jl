macro monodromy_setup(block)
    # 1. Parse input (read vars = ..., params = ... format)
    args = Dict{Symbol, Any}()
    
    # Set default values
    args[:field] = :(AcbField()) # Default field
    args[:t] = :t                # Default tracking parameter name
    
    # Parse block content
    if block.head == :block
        for line in block.args
            if line isa Expr && line.head == :(=)
                key = line.args[1]
                val = line.args[2]
                
                # Handle tuples: (x,y) -> [:x, :y]
                if val isa Expr && val.head == :tuple
                    args[key] = val.args
                elseif val isa Symbol
                     args[key] = [val] # When it's a single variable
                else
                    args[key] = val # Handle function calls like AcbField()
                end
            end
        end
    end

    # Check required inputs
    if !haskey(args, :vars)
        error("You must define 'vars = (x, y, ...)' in the setup block.")
    end

    # 2. Prepare variables
    main_vars = args[:vars]
    param_vars = haskey(args, :params) ? args[:params] : []
    t_var = args[:t]
    field_ctor = args[:field]
    
    # Internal variable η
    eta_sym = :η
    
    # 3. Generate code
    quote
        # 1. Create Field
        CCi = $(esc(field_ctor))

        # 2. Main Ring R = CC[vars, η]
        # Append η to the end of the variable list
        var_names = [$(map(string, main_vars)...), string($(QuoteNode(eta_sym)))]
        R, main_syms = polynomial_ring(CCi, var_names)
        
        # Assign user variables (x, y, z)
        $(Expr(:block, map((v, i) -> :($(esc(v)) = main_syms[$i]), main_vars, 1:length(main_vars))...))
        
        # Assign internal variable η
        $(esc(eta_sym)) = main_syms[end]

        # 3. Homotopy Ring HR = R[t]
        # Note: polynomial_ring returns the Poly object itself (not an array) when a single string is passed.
        HT, t_poly_obj = polynomial_ring(R, string($(QuoteNode(t_var))))
        $(esc(t_var)) = t_poly_obj

        # 4. Parameter Ring PR = HR[params]
        if !isempty($param_vars)
            PR, param_syms = polynomial_ring(HT, [$(map(string, param_vars)...)])
            $(Expr(:block, map((v, i) -> :($(esc(v)) = param_syms[$i]), param_vars, 1:length(param_vars))...))
            _PR = PR
        else
            _PR = HT # If no parameters, HR is the top-level ring
        end

        # Set global variables (maintain backward compatibility)
        $(esc(:_CCi)) = CCi
        $(esc(:_R))   = R
        $(esc(:_HT))  = HT
        $(esc(:_PR))  = _PR
        
        # Info messages
        println("✅ Setup Complete:")
        println("  - Variables: ", $(string.(main_vars)))
        if !isempty($param_vars)
            println("  - Parameters: ", $(string.(param_vars)))
        end
        println("  - Tracking Variable: ", $(string(t_var)))
        println("  - Internal: ", $(string(eta_sym)))
    end
end
