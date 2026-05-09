
using CertifiedHomotopyTracking

@variables x
@variables c0, c1, c2, c3, c4, c5
PREC_BITS = 256
CC = AcbField(PREC_BITS)

F = [c5*x^5+c4*x^4+c3*x^3+c2*x^2+c1*x+c0]
bp = [CC(1/2), CC(9/7), CC(9/7), CC(45/7), CC(45/56), CC(-2/3)]
x0 = [CC(.0720838,-.481488)]
v1 = vertex(bp,[x0])

vars = [x]
pars = [c0, c1, c2, c3, c4, c5]
v1 = vertex(bp,[x0])
deg = length(pars)-1

compiled_homotopy = compile_edge_homotopy(F, vars, pars)

function str_convert(
    l::Vector{Vector{Int64}}, 
    filepath::String, 
    group_name::String
)
    open(filepath, "w") do file 
        iter = 0;
        for i in l 
            write(file, "p$iter:= PermList($i);\n")
            iter = iter + 1;
        end

        write(file, "$group_name:=Group(")
        for i in 0:length(l)-1
            if i < length(l)-1
                write(file, "p$i,")
            else
                write(file, "p$i")
            end
        end
        write(file, ");")
    end
end

path_name = joinpath(@__DIR__)
filename = joinpath(path_name, "random_polynomials_5.txt")
cmd = string("Read(\"", filename, "\");")
GAP.evalstr(cmd)
@gap("G;")
@gap("StructureDescription(G);") # S5
@gap("GaloisWidth := function(G)
  local X, M, C, phi;
  if IsTrivial(G) then return 1;
  elif IsNaturalSymmetricGroup(G) or IsNaturalAlternatingGroup(G) then
    X := OrbitsDomain(G)[1];
    if Length(X) = 4 then return 3;
    else return Length(X);
    fi;
  elif IsCyclic(G) then return Maximum(Factors(Order(G)));
  elif not IsTransitive(G) then 
    return Maximum(List(Orbits(G), 
      O -> GaloisWidth(Image(ActionHomomorphism(G,O)))
    ));
  else
    X := OrbitsDomain(G)[1];
    if not IsPrimitive(G) then
      phi := ActionHomomorphism(G,Blocks(G, X),OnSets);
      return Maximum(GaloisWidth(Kernel(phi)), GaloisWidth(Image(phi)));
    elif IsSimple(G) then
      M := List(ConjugacyClassesMaximalSubgroups(G), H -> Representative(H));
      return Minimum(List(M, H -> Order(G)/Order(H)));
    else
      C := CompositionSeries(G);
      return Maximum(List([1..Length(C)-1], 
        i -> GaloisWidth(C[i]/C[i+1])
      ));
    fi;
  fi;
end;
")
@gap("GaloisWidth(G);") # 5

for n_nodes in 3:6
    result_name = "results_degree_$(deg)_$(n_nodes)_nodes_new.txt"
    result_filename = joinpath(path_name, result_name)
    path = result_filename
    gw_counts = Dict{Int, Int}()

    open(path, "w") do file
        false_count = 0;
        fail_correspondence_count = 0;
        tracking_error_count = 0;
        i = 1;
        while i <= 100
            try
                v1 = vertex(bp,[x0])
                vertices = [v1]

                for _ in 1:(n_nodes - 1) 
                    rand_u = [CC(cis(rand()*2*pi)) for _ in 1:length(pars)]
                    push!(vertices, vertex(rand_u))
                end
                edges = solve_monodromy(compiled_homotopy, vertices; max_roots=deg)


                if isempty(edges)
                    throw(ErrorException("No edges were generated."))
                end
                if length(edges[1].correspondence12) != deg
                    fail_correspondence_count += 1
                end

                dummy_name = "random_polynomials_$(deg)_$(n_nodes)nodes_new"
                save_path = joinpath(path_name, dummy_name)

                dummy_filename = joinpath(path_name, dummy_name * ".txt")
                perms = get_permutations(length(edges[1].correspondence12), edges)
                str_convert(perms, dummy_filename, "H") 
                
                cmd_dummy = string("Read(\"", dummy_filename, "\");")
                GAP.evalstr(cmd_dummy)

                A = @gap("StructureDescription(H);") 
                println(A);

                gw = -1 
                tf = string(@gap("IsomorphismGroups(G,H);"))

                if tf == "fail"
                    gw = @gap("GaloisWidth(H);")
                    gw_counts[gw] = get(gw_counts, gw, 0) + 1
                    false_count += 1
                    tf = "false"
                else
                    tf = "true"
                end
                write(file, "G#$i description: $A\n")
                write(file, "    isomorphic?: $tf\n")
                if tf == "false"
                    write(file, "    Galois width: $gw\n")
                end
                write(file, "\n")
                flush(file)
                i += 1

            catch e
                println("⚠️ Error at i=$i: $(e)")
            continue
        end
                
        end
        write(file, "number of incomplete correspondences: $fail_correspondence_count\n")
        write(file, "false count: $false_count\n")

        write(file, "Galois width counts:\n")
        for (k,v) in sort(collect(gw_counts))
            write(file, "gw = $k occurred $v times\n")
        end
    end


end


