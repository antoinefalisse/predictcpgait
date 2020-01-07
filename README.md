# predictcpgait

This repository contains code and data to generate three-dimensional muscle-driven predictive simulations of walking in children with cerebral palsy as described in: Falisse A, Pitto L, Kainz H, Hoang H, Wesseling M, Van Rossom S, Papageorgiou E, Bar-On L, Hallemans A, Desloovere K, Molenaers G, Van Campenhout A, De Groote F, and Jonkers I. Physics-based simulations to predict the differential effects of motor control and musculoskeletal deficits on gait dysfunction in cerebral palsy: a retrospective case study. Frontiers in Human Neuroscience. 2020.

Thanks for citing our work in any derived publication. Feel free to reach us for any questions: antoine.falisse@kuleuven.be | antoinefalisse@gmail.com | friedl.degroote@kuleuven.be | ilse.jonkers@kuleuven.be. This code has been developed on Windows using MATLAB2017b. There is no guarantee that it runs smooth on other platforms. Please let us know if you run into troubles.

predictcpgait consists of three main folders:

1. ParameterEstimation
    - This folder contains code to estimate personalized muscle-tendon parameters. This code is inspired from [Falisse et al. (2017)](https://ieeexplore.ieee.org/document/7748556).
2. Spasticity
    - This folder contains code to derive personalized spasticity models. This code is inspired from [Falisse et al. (2018)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0208811).
3. PredictiveSimulations
    - This folder contains code to generate predictive simulations of walking. This code is inspired from [Falisse et al. (2019a)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0217730) and [Falisse et al. (2019b)](https://royalsocietypublishing.org/doi/10.1098/rsif.2019.0402).
    
The other folders contain code and data used for the muscle-tendon parameter estimation, spasticity models, and predictive simulations:

1. CollocationScheme
    - This folder contains code for setting up direct collocation problems.
2. MuscleModel
    - This folder contains code for modeling muscle behavior and a bunch of other things.
3. OpenSimModel
    - This folder contains OpenSim model and data used in this study.
4. VariousFunctions
    - This folder contains a bunch of scripts and functions used in other code of this repository.
    