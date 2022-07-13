[![Tests and coverage](https://github.com/UofUEpi/epiworld/actions/workflows/ci.yml/badge.svg)](https://github.com/UofUEpi/epiworld/actions/workflows/ci.yml) [![C++ Documentation](https://img.shields.io/badge/Docs-Doxygen-brightgreen)](https://UofUEpi.github.io/epiworld)

# epiworld

This C++ library provides a general framework for epidemiologic simulation. The
core principle of `epiworld` is fast epidemiological prototyping for 
building complex models quickly. Here are some of its main features:

  - It only depends on the standard library (C++11 required.)
  - It is a template library.
  - It is header-only ([single file](epiworld.hpp)).
  - Models can have an arbitrary set of states.
  - Viruses and tools (e.g., vaccines, mask-wearing) can be designed to have arbitrary features.
  - Multiple tools and viruses can live in the same simulation.
  - It is *FAST*: About 30 Million person/day simulations per second (see example below).

Various examples can be found in the [examples](examples) folder.

## Hello world

Here is a simple SIR model implemented with `epiworld`. The source code
can be found [here](readme.cpp), and you can compile the code as follows:

```bash
g++ -std=c++17 -O2 readme.cpp -o readme.o
```

As you can see in [readme.cpp](readme.cpp), to use epiworld you only need
to incorporate the single header file [epiworld.hpp](epiworld.hpp):

```cpp
#include "epiworld.hpp"

using namespace epiworld;

int main()
{

    // Creating a model with three statuses:
    // - Susceptible: Status 0
    // - Infected: Status 1
    // - Recovered: Status 2
    Model<> model;
    model.add_status("Susceptible", default_update_susceptible<>);
    model.add_status("Infected", default_update_exposed<>);
    model.add_status("Recovered");

    // Desgining a virus: This virus will:
    // - Have a 90% transmission rate
    // - Have a 20% recovery rate
    // - Infected individuals become "Infected" (status 1)
    // - Recovered individuals become "Recovered" (status 2)
    // Only twenty individuals will have the virus from the beginning.
    Virus<> virus("covid 19");

    virus.set_prob_infecting(.9);
    virus.set_prob_recover(.2);
    
    virus.set_status(1, 2);

    model.add_virus_n(virus, 20);
    
    // Generating a random pop from a smallworld network
    model.population_smallworld(100000);

    // Initializing setting days and seed
    model.init(100, 123);

    // Running the model
    model.run();
    model.print();
  
}
```

And you should get something like the following:

```bash
Running the model...
_________________________________________________________________________
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| done.

________________________________________________________________________________
SIMULATION STUDY

Population size    : 100000
Days (duration)    : 100 (of 100)
Number of variants : 1
Last run elapsed t : 323.00ms
Rewiring           : off

Virus(es):
 - covid 19 (baseline prevalence: 100 seeds)

Tool(s):
 (none)

Model parameters:
 (none)

Distribution of the population at time 100:
 - (0) Susceptible :  99900 -> 34655
 - (1) Infected    :    100 -> 3159
 - (2) Recovered   :      0 -> 62186

Transition Probabilities:
 - Susceptible  0.99  0.01  0.00
 - Infected     0.00  0.86  0.14
 - Recovered    0.00  0.00  1.00
```

Which took about 0.323 seconds (~ 30 million ppl x day / second).

## Simulation Steps

The core logic of the model relies on user-defined statuses and their corresponding
update functions. In particular, the model does not have a predefined set of statuses,
e.g., susceptible, infected, recovered; it is the user who establishes them. This
provides a great deal of flexibility as models in `epiworld` can have an arbitrary
set of statuses.

Like most other ABM, `epiworld` simulates the evolution of a system in discrete steps.
Each step represents a day in the system, and changes are reflected at the beginning
of the following day. Therefore, agents can become recovered and transmit a virus
on the same day. A single step of `epiworld` features the following procedures:

1. **Status update**: Agents are updated according to the status they are at.

2. (optional) **Execute global actions**: A call of user-defined functions affecting
the system. These can make any type of change in the system.

3. (optional) **Apply rewiring algorithm**: When specified, the network is rewired
according to a user-defined function.

4. **Lock the results**: The current date is incremented in one unit and
  the changes (exposition, new infections, recoveries, etc.) are recorded
  in the database. 

5. (optional) **Mutate Variants**: When defined, variants can mutate, with the new
variants appearing the next day.

To speed up computations, `epiworld` uses by default a queuing system that decides which
agents will be active during each step and which will not. Agents are active when either
they or at least one of their neighbors has a virus active. Agents' updates are triggered
only for those who are in the queue, which in most cases accelerates the completion of
the current step.

## Agents

Agents carry two sets of important information: viruses and tools. Each agent
can have multiple instances of them, meaning that multiple viruses and tools can
coexist in a model. At each step of the simulation, an agent can face the following
changes:

- **Acquire a virus (`add_virus()`)**: Become exposed to a particular virus+host.

- **Lose a virus (`rm_virus()`)**: Removing a virus from the agent. Losing a virus
triggers a call to the virus's `postrecovery()` function, which can, for example,
result in gaining immunity to that variant. 

- **Acquire a tool (`add_tool()`)**: For example, mask-wearing, vaccines, etc.

- **Lose a tool (`rm_tool()`)**: For example, stop wearing masks, lose immunity, etc.

- **Change status (`change_status()`)**: An arbitrary change in the status of the
agent. Examples of this are moving from "exposed" to "infected," from "infected"
to "ICU," etc.

- **Become removed (`rm_agent_by_virus()`)**: An agent becomes inactive after its
condition becoming worse. In such a case, all viruses attached to the agent are
removed as well.

Any action in the model can trigger a change in its queuing system. By default,
becoming exposed makes the agent (and its neighbors) active in the queuing system.
Likewise, losing all viruses could make the agent and its neighbors inactive.

<!-- 
## Tools -->

## Contagion

Susceptible individuals can acquire a virus from any of their infected connections. The probability that susceptible individual `i` gets the virus `v` from individual `j` depends on how three things:

1. The transmissibility of the virus, <!-- $P_v \in [0,1]$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=P_v%20%5Cin%20%5B0%2C1%5D">,
2. The contagion reduction factor of `i`, <!-- $C_r \in [0,1]$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=C_r%20%5Cin%20%5B0%2C1%5D">, and
3. The host's transmission reduction factor, <!-- $T_r \in [0,1]$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=T_r%20%5Cin%20%5B0%2C1%5D">.

The last two are computed from <!-- $i$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=i"> and <!-- $j$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=j">'s tools. Ultimately, the probability of <!-- $i$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=i"> getting virus $v$ from <!-- $j$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=j"> equals:

<!-- $$
P(Virus v) = P_v * (1 - C_r) *(1 - T_r) 
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=P(Virus%20v)%20%3D%20P_v%20*%20(1%20-%20C_r)%20*(1%20-%20T_r)%20"></div>

Nonetheless, the default behavior of the simulation model is to assume that individuals can acquire only one disease at a time, if any. This way, the actual probability is:

<!-- $$
P(Virus v|\mbox{at most one virus}) = P_{ivj}
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=P(Virus%20v%7C%5Cmbox%7Bat%20most%20one%20virus%7D)%20%3D%20P_%7Bivj%7D"></div>

The latter is calculated using Bayes' rule

<!-- $$
\begin{align*}
P_{ivj} & = P(at most one virus|Virus v) * P_v / P(at most one virus) \\
        & = P(Only Virus v)/P_v * P(Virus v) / P(at most one virus) \\
        & = P(Only Virus v)/P(at most one virus)
\end{align*}
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0AP_%7Bivj%7D%20%26%20%3D%20P(at%20most%20one%20virus%7CVirus%20v)%20*%20P_v%20%2F%20P(at%20most%20one%20virus)%20%5C%5C%0A%20%20%20%20%20%20%20%20%26%20%3D%20P(Only%20Virus%20v)%2FP_v%20*%20P(Virus%20v)%20%2F%20P(at%20most%20one%20virus)%20%5C%5C%0A%20%20%20%20%20%20%20%20%26%20%3D%20P(Only%20Virus%20v)%2FP(at%20most%20one%20virus)%0A%5Cend%7Balign*%7D"></div>


Where

<!-- $$
\begin{align*}
P(Only Virus V)      & = P_v * \prod_{m!=V} (1 - P_m) \\
P(at most one virus) & = P(None) + \sum_{k \in viruses} P_k * \prod_{m != k} (1 - P_m) \\
P(None)              & = \prod_{k \in Viruses} (1 - P_k)
\end{align*}
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0AP(Only%20Virus%20V)%20%20%20%20%20%20%26%20%3D%20P_v%20*%20%5Cprod_%7Bm!%3DV%7D%20(1%20-%20P_m)%20%5C%5C%0AP(at%20most%20one%20virus)%20%26%20%3D%20P(None)%20%2B%20%5Csum_%7Bk%20%5Cin%20viruses%7D%20P_k%20*%20%5Cprod_%7Bm%20!%3D%20k%7D%20(1%20-%20P_m)%20%5C%5C%0AP(None)%20%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%3D%20%5Cprod_%7Bk%20%5Cin%20Viruses%7D%20(1%20-%20P_k)%0A%5Cend%7Balign*%7D"></div>

This way, viruses with higher transmissibility will be more likely to be acquired when competing with other variants.
