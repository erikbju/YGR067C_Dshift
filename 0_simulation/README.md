# Simulation

## Setting up environment

We recommend that in order to run the simulation experiments that
you use the Docker configuration provided. To do so, simply run
the following command from this directory.

```sh
docker compose run --build --rm simulations
```

## Running the experiment

ALl the settings for the experiment are passed to the simulation
script via a `.yaml` file. To run the experiment with the same
parameters as described in the manuscript, you can use the file
provided.

```sh
python simulation_experiment.py ygr067c-experiment-settings.yaml
```

