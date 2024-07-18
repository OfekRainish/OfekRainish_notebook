
## Statistical Analysis of Chordata Presence in Mediterranean Sea Sites

In this statistical analysis, I aim to investigate the occurrence percentage of Chordata across various depths and seasons at two Middle Eastern sites: Achziv and Sdot Yam. The primary objective is to determine if there exists a statistical relationship between these variables and the presence of Chordata.

#### Methods

I will employ three statistical tools:

1. **Bartlett Test**: This will assess whether the variances of the Chordata variable are homogeneous across the different sites.
   
2. **Kruskal-Wallis Test**: This test will evaluate if the distribution of Chordata significantly differs among the different sites.
   
3. **Type III ANOVA Tests**: These tests will perform an analysis of variance on the Chordata variable with respect to site, depth, season, and their interactions.


### Results
1. When I conducted a bartlett test to examine whether the variances of the Chordata variable are equal across the different sites. I got a p value of **2.75e-12**, which indicates that the variances of chordata are significantly different across the different sites.

2. When conducting a Kruskal-Wallis test to assess whether the distribution of Chordata varies significantly across different sites, a p-value of **6.981e-07** was obtained. This suggests a statistically significant difference in the distribution of Chordata among the two sites.


3.  Anova Table

| Term                         | Df  | Df.res | F value    | Pr(>F)      | Significance |
|------------------------------|-----|--------|------------|-------------|--------------|
| site                         | 1   | 328    | 47.055287  | 3.4464e-11  | ***          |
| as.factor(depth)             | 2   | 328    | 68.389343  | < 2.22e-16  | ***          |
| season                       | 1   | 328    | 0.038526   | 0.84451     |              |
| site:as.factor(depth)        | 2   | 328    | 16.051637  | 2.2350e-07  | ***          |
| site:season                  | 1   | 328    | 0.019229   | 0.88980     |              |
| as.factor(depth):season      | 2   | 328    | 12.861264  | 4.1940e-06  | ***          |
| site:as.factor(depth):season | 2   | 328    | 17.176740  | 8.0462e-08  | ***          |

 Conclusions:


   -  The site has a highly significant effect on Chordata.

   -  Depth has an extremely significant effect on Chordata.

   -  There is no significant effect of the season on Chordata.

   -  The interaction between site and depth is highly significant, meaning **the effect of the site on Chordata depends on the depth**.

   - There is no significant interaction effect between site and season on Chordata.

   -  The interaction between depth and season is highly significant, meaning **the effect of depth on Chordata depends on the season**.

   -  The three-way interaction between site, depth, and season is highly significant. This means the combined effect of site and depth on Chordata varies depending on the season.

4. 

![alt text](../images/data%20anlysisn%20class.png)


   It can be seen from the graph that chordata were not found at all at a depth of 10 meters in any of the sites.
Another interesting thing is that at the Akziv site, Chordata were only found at a depth of 45 meters regardless of the season.

## Correlation between Porifera and Chordata
Here I wanted to check if there is a correlation between these two groups.

![alt text](../images/cor%20plot%201.png)

The correlation coefficient is **0.3116275**. It means that there is a weak to moderate tendency for Porifera and Chordata to increase together, according tp the Photosurvey dataset. However, this connection is not very strong, so there is still a lot of variation that this relationship doesn't explain.

### Supplements


[raw data](../exel%20files/Photosurvey_metadata%20class.csv)