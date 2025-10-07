#!/bin/bash

COMMANDS=(
    "gcc -Wall -Iinclude -Ilib -g -O0 -o bin/main_no_learn_simulation_cuda src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c src/neuron_models/GPU_lif_neuron.cu -lm -DNO_LEARN -DOPENMP -DCUDA"
)

#COMMANDS=(
#    "gcc -Wall -Iinclude -Ilib -o bin/main_parallel_no_learn_simulation src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm -DOPENMP -DNO_LEARN -fopenmp"
#    "gcc -Wall -Iinclude -Ilib -o bin/main_parallel_no_learn_simulation_reorder src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm -DOPENMP -DREORDER -DNO_LEARN -fopenmp"
#    "gcc -Wall -Iinclude -Ilib -o bin/main_no_learn_simulation src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm -DNO_LEARN -DOPENMP"
#    "gcc -Wall -Iinclude -Ilib -o bin/main_no_learn_simulation_reorder src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm -DREORDER -DNO_LEARN -DOPENMP"
#    "gcc -Wall -Iinclude -Ilib -g -O0 -o bin/main_parallel_no_learn_simulation_reorder_g src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm -DREORDER -DNO_LEARN -DOPENMP"
#    "gcc -Wall -Iinclude -Ilib -g -O0 -o bin/main_parallel_no_learn_simulation_g src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm -DNO_LEARN -DOPENMP"
#    "gcc -Wall -Iinclude -Ilib -g -O0 -o bin/main_no_learn_simulation_cuda src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c src/neuron_models/GPU_lif_neuron.cu -lm -DNO_LEARN -DOPENMP -DCUDA"
#)




#    


#COMMANDS=(
#    "gcc -Iinclude -Ilib -o bin/main_parallel_learn_simulation src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DOPENMP"
#    "gcc -Iinclude -Ilib -o bin/main_parallel_no_learn_simulation src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DOPENMP -DNO_LEARN"
#    "gcc -Iinclude -Ilib -o bin/main_parallel_simulation_debug src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DOPENMP -DDEBUG"
#    "gcc -Iinclude -Ilib -o bin/main_simulation src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1"
#    "gcc -Iinclude -Ilib -o bin/main_simulation_no_learn src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DNO_LEARN"
#    "gcc -Iinclude -Ilib -o bin/main_simulation_debug src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DDEBUG"
#    "gcc -Iinclude -Ilib -g -o bin/main_parallel_learn_simulation_gdb src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DOPENMP"
#    "gcc -Iinclude -Ilib -o bin/main_simulation src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DOPENMP -DNOLEARN"
#    "gcc -Iinclude -Ilib -o bin/main_simulation_reorder src/main.c src/snn_library.c src/load_data.c src/helpers.c src/neuron_models/lif_neuron.c src/training_rules/stdp.c src/simulations/simulations.c -lm  -DINPUT_SYNAPSES=1 -DINPUT_NEURON_BEHAVIOUR=2 -DINPUT_WEIGHTS=1 -DINPUT_DELAYS=2 -DINPUT_TRAINING_ZONES=1 -DOPENMP -DNOLEARN -DREORDER"
#
#    "gcc -Iinclude -Ilib -o bin/network_generator src/network_generator/network_generator_main.c src/network_generator/network_generator.c -lm"
#
#    "gcc -Iinclude -Ilib -o bin/image_encoder src/encoders/image_encoder_main.c src/encoders/image_encoders.c -lm"
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