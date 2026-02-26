#!/bin/bash


COMMANDS=(
"gcc -O3 -ffast-math -fopenmp -funroll-loops -fprefetch-loop-arrays -flto -mtune=native -Iinclude -Iinclude/pultsito -Ilib -o bin/main_simulation_opt src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c -lm"
"gcc -O3 -march=skylake-avx512 -ffast-math -funroll-loops -fprefetch-loop-arrays -flto -mtune=native -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation_opt src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c -lm -DAVX512"
"gcc -Ilib -o bin/snngenerator src/network_generator/snn_generator.c -lm"
)

#COMMANDS=(
#"gcc -O3 -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -g -o bin/main_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm"
#"gcc -O3 -march=skylake-avx512 -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DAVX512"
#"gcc -O3 -march=skylake-avx512 -ffast-math -funroll-loops -fprefetch-loop-arrays -flto -mtune=native -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation_opt src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DAVX512"
#"gcc -O3 -ffast-math -fopenmp -funroll-loops -fprefetch-loop-arrays -flto -mtune=native -Iinclude -Iinclude/pultsito -Ilib -o bin/main_simulation_opt src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm"
#)

echo "== Compiling... =="

# ejecutar cada comando
for CMD in "${COMMANDS[@]}"; do
    echo "  Generating executable..."
    # ejecutar el comando
    eval $CMD
done

echo ""
echo "== Code compiled! =="