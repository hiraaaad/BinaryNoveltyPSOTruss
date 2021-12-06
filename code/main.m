function main(benchmark)

switch benchmark
    case {10, 25, 52}
        output = exact_enumeration(benchmark);
    case {15, 224, 68, 47, 72, 200}
        output = ndbpso(benchmark);
end

output