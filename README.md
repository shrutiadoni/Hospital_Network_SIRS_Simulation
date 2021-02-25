# Hospital Network SIRS Simulation
SIRS- (Susceptible, Infected, Removed and Susceptible again)
Dataset: http://www.sociopatterns.org/datasets/hospital-ward-dynamic-contact-network/


Social interactions, besides providing experiences, satisfactions and joy between people can also be the easiest way for diseases to transmit among others. Since these interactions are somehow unavoidable, the risk of infections of some diseases is indeed, unavoidable as well.

We will conduct network statistical analysis using Gephi and R, as well as epidemic modeling using two different models in R. We used the same edges information and created a graph to induce SIRS infection/epidemic. Here we created a function to track the status of each node and four vectors- Susceptible, Infected, and removed.

Infection is induced based on the node’s links i.e. interaction with other nodes.

Experiment:
- For experimentation purpose, we considered remove after as 2 and susceptible again as 3. We chose these as a hospital is easy place for recovery and quarantine and it enabled proper observation of all statuses of nodes.
- This means once a node is infected, after time +2, it gets recovered/removed
- After removal, once t=t+3, it becomes susceptible again.
