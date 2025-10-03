#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


// Function to generate random spike trains
int random_spike_train_generator(int *spike_train, int time_steps, int prob){
    
    int p, t, n_spks = 0;

    // generate spikes 
    for(t=0; t<time_steps; t++){
        
        p = rand() % 100;

        // generate spike on time t if p > prob
        if(p > prob){
            
            spike_train[n_spks] = t;
            n_spks ++;
        }
    }

    // return number of spikes
    return n_spks;
}


int main(int argc, char *argv[]) {

    int i, j;
    int n_input_spike_trains, time_steps, prob;
    int *spk_train, n_spikes;
    FILE *f;

    
    // randomize
    srand(time(NULL)); 

    // read the number of spike trains to generate
    n_input_spike_trains = strtoul(argv[1], NULL, 10);

    // read the time steps for each spike train
    time_steps = strtoul(argv[2], NULL, 10);

    // read the probability to generate a spike in a time step
    prob = strtoul(argv[3],NULL, 10);

    // open file
    f = fopen(argv[4], "w"); //f = fopen("input_spikes.txt", "w");
    if(f == NULL){
        printf(" > Error opening the file!\n");
        exit(1);
    }

    // allocate memory for spike trains
    spk_train = (int *)malloc(time_steps * sizeof(int));


    // generate spike trains 
    for(i = 0; i<n_input_spike_trains; i++){
        
        // generate random spike train
        n_spikes = random_spike_train_generator(spk_train, time_steps, prob);

        // store spike train
        fprintf(f, "%d ", n_spikes);

        for(j = 0; j<n_spikes; j++){
            fprintf(f, "%d ", spk_train[j]);
        }
        fprintf(f, "\n");
    }

    // close file
    fclose(f);
}