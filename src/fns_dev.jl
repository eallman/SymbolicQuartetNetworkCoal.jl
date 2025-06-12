# Working on improvements

# List of functions in this file
#
# function network_expectedCF_formulas_v2(network::HybridNetwork; 
#         showprogressbar=false, 
#         inheritancecorrelation=0, 
#         filename="symbolicQNC_HFO_out"::AbstractString,     
#         symbolic=false::Bool,
#         savecsv=false::Bool,
#         macaulay=false::Bool,
#         matlab=false::Bool,
#         multigraded=false::Bool,
#         singular=false::Bool,
#         dpoints::Integer=7   # ESA
#         )
# function printCFs(qData::Vector{PN.QuartetT{MVector{3, Float64}}},
# function makeEdgeLabel_v2(net::PN.HybridNetwork; showAllEdgeLabels::Bool=false)

#=
tauSymbol = "t_"

# Dictionary for edge lengths (τ)
translationBLs = Dict()
numEdges = length(net.edge)
allNetEdges = collect(1:numEdges)

termEdges = [e.number for e in net.edge if PhyloNetworks.getchild(e).leaf]
intNetEdges = setdiff(allNetEdges, termedgenum)

for e in intNetEdges
    translationBLs[net.edge[e].length] = "$tauSymbol{$e}"
end

translationBLs = Dict(value => key for (key, value) in translationBLs)
translationBLs

##
for ln = 1:11
    net_newick = readline("output.log")
    println(net_newick)
end
##

# alllines = readlines("output.log")

=#

# search on symCF and ESA
"""
   network_expectedCF_formulas_v2(network::HybridNetwork; 
                                showprogressbar=true, 
                                inheritancecorrelation=0, 
                                filename="symbolicQNC-HFO-out"::AbstractString, 
                                symbolic=false::Bool, 
                                savecsv=false::Bool, 
                                macaulay=false::Bool, 
                                matlab=false::Bool, 
                                multigraded=false::Bool,
                                singular=false::Bool)

Compute expected concordance factors (CFs) for quartets in a phylogenetic network, optionally generating symbolic formulas.

This function calculates the expected concordance factors (CFs) for all possible quartets of taxa in a `HybridNetwork`, 
taking into account coalescent processes and inheritance correlations at hybrid nodes. It can operate in numerical mode 
(default) or symbolic mode (when `symbolic=true`), where it expresses CFs as formulas involving branch lengths and 
inheritance parameters (γ). Results are logged to a file, and optional outputs can be saved as CSV, Macaulay2, MATLAB, 
or multigraded implicitization files for further analysis.

# Arguments
- `network::HybridNetwork`: A phylogenetic network from the PhyloNetworks package, with edge lengths in coalescent units and γ values for hybrid edges.
- `showprogressbar::Bool=true`: If true, displays a progress bar for quartet calculations.
- `inheritancecorrelation::Number=0`: Correlation between inheritance probabilities at hybrid nodes (must be between 0 and 1).
- `filename::AbstractString="symbolicQNC_HFO_out"`: Base name for output files (e.g., log, CSV, Macaulay2, MATLAB).
- `symbolic::Bool=false`: If true, computes CFs as symbolic expressions; requires all edge parameters to be defined or assigns random values.
- `savecsv::Bool=false`: If true, saves CFs to a CSV file named `<filename>.csv`.
- `macaulay::Bool=false`: If true, generates a Macaulay2 script (`<filename>.m2`) for symbolic analysis; requires `symbolic=true`.
- `matlab::Bool=false`: If true, generates a MATLAB script (`<filename>.m`) for symbolic analysis; requires `symbolic=true`.
- `multigraded::Bool=false`: If true, generates a Macaulay2 multigraded implicitization script (`<filename>.im.m2.txt`); requires `symbolic=true`.
- `singular::Bool=false`: If true, generates a Singular script (`<filename>.sing`); requires `symbolic=true`,
- `dpoints::Integer=7`: Number of digits of accuracy for numerical CFs .

# Returns
- `quartet::Vector{PhyloNetworks.QuartetT}`: Array of QuartetT objects containing quartet indices, taxa, and CFs (numerical or symbolic).
- `taxa::Vector{String}`: Sorted list of taxon names from the network.

# Throws
- `ErrorException`: If the root is a leaf, edge lengths are missing/negative, γ values are missing/negative, or inheritance correlation is invalid.
- `ErrorException`: If `macaulay`, `matlab`, or `singular` is true but `symbolic` is false.

# Side Effects
- Writes a log file (`<filename>.log`) with topology, parameters, and CFs.
- Optionally writes CSV, Macaulay2, MATLAB, multigraded, or Singular files based on flags.
"""
function network_expectedCF_formulas_v2(network::HybridNetwork; 
                            showprogressbar=false, 
                            inheritancecorrelation=0, 
                            filename="symbolicQNC_HFO_out"::AbstractString,     
                            symbolic=false::Bool,
                            savecsv=false::Bool,
                            macaulay=false::Bool,
                            matlab=false::Bool,
                            multigraded=false::Bool,
                            singular=false::Bool,
                            dpoints::Integer=7   # ESA
                            )
    
    #----------filename and logging----------#
    log=string(filename,".log")
    logfile=open(log,"w")
    str="SymbolicQuartetNetworkCoal.jl log\n"
    str*="Timestamp: $(Dates.now())\n"
    str*="------------------------\n"
    write(logfile,str)
    flush(logfile)
    #---------forcing things to work---------#
    str="General setting: \n"
    str*="Symbolic option: "
    if(symbolic) 
        net=deepcopy(network)        
        #check if params all exist otherwise, assign them using readTopologyrand
        params=true
        for e in net.edge
            if e.length<0
                params=false
            end
        end

        if !(params)
            @warn("Input topology is missing some parameters. Assigning arbitrary values using [readTopologyrand]")
            net=readTopologyrand(net)
        else
            net=network
        end

        str*="on\n"
    else
        net=network 
        str*="off\n"
    end
    str*="Store output in .csv file: "
    str*=(savecsv ?  "on\n" : "off\n")
    if(matlab) symbolic || error("symbolic must be set to true.") end
    str*="Write Matlab file: "
    str*=(matlab ? "on\n" : "off\n")
    
    if(macaulay) symbolic || error("symbolic must be set to true.") end
    str*="Write Macaulay2 file: "
    str*=(macaulay ? "on\n" : "off\n")

    if(singular) symbolic || error("symbolic must be set to true.") end
    str*="Write Singular file: "
    str*=(macaulay ? "on\n" : "off\n")

    write(logfile,str)
    flush(logfile)

    #-----------topologies and dictionary-----------#
    # ESA ESA
    dict=symCF.parameterDictionary(net,inheritancecorrelation)
    str="------------------------\n"
    str*="Topology:\n"
    str*="$(writeTopology(net,digits=dpoints))\n"
    if(symbolic) str*="$(symCF.gettingSymbolicTopology(net,dict))\n" end
    write(logfile,str)
    flush(logfile)

    # ESA ESA
    # ----------write dictionary to file-----------#
    if (symbolic)
        dict_fname=string(filename,".dict")
        dfile=open(dict_fname,"w")
        for e in net.edge
            str *= "$(dict[e.length])\t\t$(e.length)\n"
            if (e.hybrid)
                str *= "$(dict[e.gamma])\t\t$(e.gamma)\n"
            end
        end
    end

    # ESA ESA end

    #-----------setting up desk-----------#
    #need to fix when some parameter values are identical
    df = DataFrame(Split=String[], CF=String[]) 
    str="------------------------\n"
    str*="Parameters:\n"
    str*="Parameter\t\tValue\n"
    for e in net.edge
        str*="$(dict[e.length])\t\t$(e.length)\n"
        if(e.hybrid) str*="$(dict[e.gamma])\t\t$(e.gamma)\n" end
    end
    str*="$(dict[inheritancecorrelation])\t\t$(inheritancecorrelation)\n"
    if symbolic
        write(logfile,str)
        flush(logfile)    
    else
        str=""
        write(logfile,str)
        flush(logfile)    
    end
    
    #--------(almost) original stuff--------#
    net.node[net.root].leaf && error("The root can't be a leaf.")
    PN.check_nonmissing_nonnegative_edgelengths(net,
        "Edge lengths are needed in coalescent units to calcualte expected CFs.")
    all(e.gamma >= 0.0 for e in net.edge) || error("some γ's are missing for hybrid edges: can't calculate expected CFs.")
    inheritancecorrelation >= 0 || error("the inheritance correlation should be non-negative")
    inheritancecorrelation <= 1 || error("the inheritance correlation should be <= 1")
    taxa = sort!(tipLabels(net))
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = PN.nchoose1234(ntax) # matrix to rank 4-taxon sets
    qtype = MVector{3,Float64} # 3 floats: CF12_34, CF13_24, CF14_23; initialized at 0.0
    numq = nCk[ntax+1,4]
    quartet = Vector{PN.QuartetT{qtype}}(undef, numq)
    ts = [1,2,3,4]
    for qi in 1:numq
        quartet[qi] = PN.QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.))
        # next: find the 4-taxon set with the next rank,
        #       faster than using the direct mapping function
        ind = findfirst(x -> x>1, diff(ts))
        if ind === nothing ind = 4; end
        ts[ind] += 1
        for j in 1:(ind-1)
            ts[j] = j
        end
    end
    if showprogressbar
        nstars = (numq < 50 ? numq : 50)
        nquarnets_perstar = (numq/nstars)
        println("Calculation quartet CFs for $numq quartets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end
    for qi in 1:numq
        symCF.network_expectedCF!(quartet[qi], net, taxa, taxonnumber, inheritancecorrelation, df, symbolic, dict)
        if showprogressbar && qi >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * nquarnets_perstar))
        end
    end
    showprogressbar && print("\n")

    #--------output--------#
    str="------------------------\n"
    str*="Concordance factor:\n"
    str*="Quartet\t\tFormula\n"
    for i in 1:nrow(df)
        str*="$(df[i,1])\t\t$(df[i,2])\n"
    end
    if(savecsv) CSV.write("$filename.csv",df,header=false) end
    write(logfile,str)
    flush(logfile)
    
    #macaulay output
    numCFs=size(df)[1]
    if(symbolic) 
        dataframe=deepcopy(df)
        params=symCF.gettingSymbolicInput(net, dataframe, inheritancecorrelation) 
    end
    if(macaulay)
        open("$filename.m2", "w") do file
        str="R = QQ["
        for par in params str=str*par*"," end
        str=str*"C_1..C_$numCFs]\n"
        str=str*"I = ideal(\n"
        i=1
        while i<numCFs
            str=str*"$(dataframe[i,2])-C_$i,\n"
            i+=1
        end
        str=str*"$(dataframe[numCFs,2])-C_$numCFs);\n"
        str=str*"G = eliminate (I, {"
        numparams=length(params)   
        for par in 1:(length(params)-1) str=str*params[par]*"," end
        str=str*"$(params[numparams])})\n"
        str=str*"S = QQ[C_1..C_$numCFs]\nJ = sub(G,S)\ndim J"
        write(file, str)
        end
    end

    #matlab output
    if(matlab)
        open("$filename.m", "w") do file 
            str="% Declare variables\n"
            str=str*"syms "
            for par in params str=str*par*" " end
            for i in 1:numCFs str=str*"C_$i " end
            str=str*"\n\n% matrix of generating polynomials\n"
            str=str*"F=["
            i=1
            while i<numCFs
                str=str*"$(dataframe[i,2])-C_$i,\n"
                i+=1
            end
            str=str*"$(dataframe[numCFs,2])-C_$numCFs];\n"
            str=str*"\n% matrix of generating polynomials\n"
            str=str*"\n% Array of all variables\n"
            str=str*"V=["
            for par in params str=str*par*" " end
            for i in 1:numCFs-1 str=str*"C_$i " end
            str=str*"C_$numCFs]\n"
            str=str*"\n% Compute dimension\nCoalDim(F,V)"
        write(file, str)
        end
    end

    if(multigraded)
        open("$filename.im.m2.txt", "w") do file
        str="needsPackage \"MultigradedImplicitization\"\n"
        str*="R = QQ["
        for par in params str=str*par*"," end
        str*="T]\n"
        str*="S = QQ["
        str=str*"C_0..C_$numCFs]\n"
        str=str*"im = {\n"
        i=1
        while i<numCFs
            str=str*"$(dataframe[i,2]),\n"
            i+=1
        end
        str=str*"$(dataframe[numCFs,2])};\n"

        str=str*"im = {T} | apply(im, f -> f * T);\n"
        str=str*"phi = map(R, S, im)\n"
        str=str*"d=1\n"
        str=str*"I = time componentsOfKernel(d, phi)\n"
        str=str*"L = flatten values I"

        write(file, str)
        end
    end

    if(singular)
        open("$filename.sing", "w") do file
        str="ring R = 0, ("
        for par in params str=str*par*"," end
        for i in 1:numCFs str=str*"C$i," end
        str=chop(str)
        str=str*"), dp;\n"
        str=str*"ideal I = \n"
        i=1
        while i<numCFs
            str=str*"$(dataframe[i,2])-C$i,\n"
            i+=1
        end
        str=str*"$(dataframe[numCFs,2])-C$numCFs;\n"

        str=str*"ideal G = eliminate (I, "
        for par in params str=str*par*"*" end
        str=chop(str)
        str=str*");\n"

        str=str*"ring S = 0, ("
        for i in 1:numCFs str=str*"C$i," end
        str=chop(str)
        str=str*"), dp;\n"
        str=str*"ideal J = imap(R, G);\n"
        str=str*"int dimJ = dim(J);\n"
        str=str*"dimJ;"

        write(file, str)
    end
end

    return quartet, taxa, dict
end

"""
    printCFs(qData::Vector{PN.QuartetT{MVector{3, Float64}}},
             tData::Vector{String};
             numDigits::Integer=4)

# Pretty prints quartet CFs using numDigits
    
"""
function printCFs(qData::Vector{PN.QuartetT{MVector{3, Float64}}},
        tData::Vector{String}; 
        numDigits::Integer=4)

    for qt in qData
        print(tData[qt.taxonnumber])
        print(" ")
        println(map(x -> round(x,digits=numDigits),qt.data))
    end
    return
end

"""
    makeEdgeLabel_v2(net; showTerminalEdgeLabels=false)

Generates a dataframe mapping edge numbers to their symbolic labels.

## Description
This function creates labels for the edges of a `HybridNetwork` in the format `"t_{e}"`, where `e` is the edge number.  
By default, labels are only assigned to **non-terminal edges** (i.e., edges that do not end at leaf nodes).  
The output dataframe is used as input for PhyloPlots' option `edgelabel=``.
Setting `showTerminalEdgeLabels=true` includes labels for terminal edges as well.

## Arguments
- `net`: A `HybridNetwork` object.
- `showAllEdgeLabels`: A boolean flag (default = `false`).  
   - `false`: Excludes terminal edges, hybrid edges with one leaf descendant, 
                and terminal edges with parent the root.  
   - `true`: Includes all edges.  

## Returns
- A `DataFrame` with columns:
  - `number`: Edge numbers.
  - `label`: Corresponding symbolic labels (`"t_{e}"`).
"""
function makeEdgeLabel_v2(net::PN.HybridNetwork; showAllEdgeLabels::Bool=false)

    # ESA update to avoid branch parameters which do not appear in symbolic CFs

    # get internal edge numbers unless want all edges labeled
    edge_numbers_to_include = [e.number for e in net.edge if !PhyloNetworks.getchild(e).leaf || showAllEdgeLabels]

    if !showAllEdgeLabels

        # possibly remove some edges from list if not needed in symbolic parameterization
        hybridEdgeInds = findall(n -> n.hybrid, net.edge)
        hybridEdgeNumbers = map(e -> e.number, net.edge[hybridEdgeInds])
        hybridNodeInds = findall(n -> n.hybrid, net.node)

        # ESA
        println(hybridEdgeInds)
        println(hybridEdgeNumbers)
        println(hybridNodeInds)

        edge_numbers_to_remove = []

        # first check if root is parent of leaf
        rootHasLeafChild = false
        for e in net.node[net.root].edge
            if getchild(e).leaf
                rootHasLeafChild = true
                edge_numbers_to_remove = [ee.number for ee in net.node[net.root].edge]
                break
            end
        end

        # now check if hyrid nodes have only one descendant
        hybridNodesWithOneLeafDescendant = [n for n in net.node[hybridNodeInds] if getchild(n).leaf]
        # ESA
        # PP.plot(net, shownodenumber=true)
        # println(" ***** ")
        # println(hybridNodesWithOneLeafDescendant)
        # println(" ***** ")
        foreach(x -> push!(edge_numbers_to_remove, getparentedge(x).number), hybridNodesWithOneLeafDescendant)
        foreach(x -> push!(edge_numbers_to_remove, getparentedgeminor(x).number), hybridNodesWithOneLeafDescendant)

        if !isempty(edge_numbers_to_remove)
            setdiff!(edge_numbers_to_include, edge_numbers_to_remove)
        end
    end

    println(edge_numbers_to_include)

    df = DataFrame(
        number=[num for num in edge_numbers_to_include],
        label=["t_{$num}" for num in edge_numbers_to_include]
    )

    return df
end



# ##
# R"par(mfrow=c(1,2))"
# ##
# plot(ntwk,showedgenumber=true);
# plot(ntwk,shownodenumber=true);
# ##