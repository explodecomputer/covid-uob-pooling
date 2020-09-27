# Pooling strategies for detecting Covid in University of Bristol students

The University of Bristol UoB plans to commence the 2020 academic year with a strategy of placing students in living circles, to minimise the spread of SARS-CoV-2 across the student body. To monitor outbreaks it would be costly to simply test everybody on a regular basis, and so the question of whether pooling samples could effectively reduce costs and still adequately detect outbreaks is an open question.

These simulations seek to explore pooling strategies that aim to 

- maximise sensitivity (maximise true positive detection rate)
- maximise specificity (minimise true negative detection rate)
- minimise cost
- minimise reagent use

## Infections in the population

The simulations will generate a population the size of the student body, and place them into hierarchies that represent location and living circle. For a given prevalence, sample individuals to be positive based on clustering within these hierarchies. The less clustered the positive cases are would represent less adherence to social distancing in the student population.

## Assay sensitivity

Sensitivity is assumed to be affected by pooling. The assay for detecting SARS-CoV-2 will perform worse if the dilution of the sample is greater, so a larger pool with one case will have lower sensitivity than a smaller pool with one case.

Assume 98% sensitivity and 100% specificity for a single sample (not pooled / diluted). By pooling samples, the dilution is assumed to be the proportion positive divided by the proportion in the pool. In order to simulate this appropriately we need a function that will relate sensitivity to dilution fraction.

## Variables within the simulation

- Baseline - test everybody
- Larger vs smaller pools
    + Larger pools will have the same amount of time
    + Larger pools will cost proportionately less
    + Larger pools will have lower sensitivity
    + Larger pools will use less reagent
- Detecting infected living circles vs infected individuals
    + Detecting just living circles will have lower cost
    + Detecting just living circles will have lower specificity (not everyone in the living circle will be infected, track and trace)
    + Detecting just living circles will have lower reagent use

## Outbreak

Define outbreak as 2 people in a location being infected. This is detected either from

- individual tests returning 2 people in a location
- one location pool turning up positive
- pool + followup tests returning 2 people in a location



## Other thoughts

If a pool is detected positive, is it possible to avoid having to test everyone in the pool to identify the infected individuals?

Group testing theory - when prevalence is low then pooling can be more efficient, basically just look for positive results in a pool and then follow up everyone in the pool. To maximise the efficiency of the pooling need to have a sense of the prevalence. 

This paper uses tags along with pooling:
https://advances.sciencemag.org/content/6/37/eabc5961

This is an R package for group testing:
https://journal.r-project.org/archive/2010/RJ-2010-016/RJ-2010-016.pdf


https://linkinghub.elsevier.com/retrieve/pii/S1198743X20303499

https://theconversation.com/group-testing-for-coronavirus-called-pooled-testing-could-be-the-fastest-and-cheapest-way-to-increase-screening-nationwide-141579

