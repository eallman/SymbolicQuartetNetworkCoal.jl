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


function parameterDictionary_XX(net, inheritancecorrelation; tauSymbol::String="t_", gammaSymbol::String="r_")
    dict = Dict()
    # Dictionary for edge lengths (τ)
    numEdges = length(net.edge)
    allNetEdges = collect(1:numEdges)
    termedgenum = [e.number for e in net.edge if PhyloNetworks.getchild(e).leaf]
    maxMergedEdges = net.numTaxa + 2*(net.numHybrids) #max_edges_from_leaf_to_leaf(net)-2
    for mergeRange in 1:maxMergedEdges
        if mergeRange == 1
            for e in allNetEdges
                dict[net.edge[e].length] = "$tauSymbol{$e}"
            end
        else
            edgeCombinations = collect(combinations(allNetEdges, mergeRange))
            for edgeCombo in edgeCombinations
                if isempty(intersect(edgeCombo, termedgenum))
                    length, symbolicName = mergedEdgeLengthandSymbolicName(net, edgeCombo)
                    dict[length] = symbolicName
                end
            end
        end
    end
    # Dictionary for inheritance probabilities (γ)
    hybridNodeNumbers = [n.number for n in net.node if n.hybrid]
    for j in 1:net.numHybrids
        hybNode = hybridNodeNumbers[j]
        visitCount = 1
        for e in net.edge
            if PhyloNetworks.getchild(e).number == hybNode
                e.gamma = round(e.gamma, digits=dpoints)
                if visitCount == 1
                    dict[e.gamma] = "$gammaSymbol{$j}"
                    visitCount += 1
                elseif visitCount == 2
                    dict[e.gamma] = "(1-$gammaSymbol{$j})"
                else
                    error("Hybrid node $(hybNode) has more than 2 incoming edges.")
                end
            end
        end
    end
    # Inheritance correlation (ρ)
    inheritancecorrelation = round(inheritancecorrelation, digits=dpoints)
    dict[inheritancecorrelation] = "&rho"
    dict[round(1 - inheritancecorrelation, digits=dpoints)] = "1-&rho"
    return dict
end