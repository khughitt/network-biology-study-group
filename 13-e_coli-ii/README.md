Many Microbe Microarrays Database (M3D) Data Preparation
========================================================

<a href='mailto:khughitt@umd.edu'>Keith Hughitt</a> (<time>2014-07-17</time>)

[view source](README.rmd)

Overview
========

Below we will load in some E. coli transcriptome data from the [Many Microbe Microarrays Database (M<sup>3D</sup>)](http://www.m3d.mssm.edu), and explore some basic properties of the data.

Hopefully, the code below will provide a good starting point for future R-based network analyses of the data.

The data included in this repo and used below comes M<sup>3D</sup> build 6, accessed on July 16, 2014. Currently, the gene only normalized version of the **E. coli** microarray data is used. In the future it may be worth analyzing the raw data and including probes for intergenic regions.

Methods
=======

Load libraries
--------------

``` {.r}
library(cbcbSEQ)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(gplots)
library(ggplot2)
library(limma)
```

Load data
---------

### Load feature data

To begin, let's load in some of the metadata which describes the sources and experimental conditions included in the microarray datasets.

``` {.r}
# load feature data
features_long = tbl_df(read.delim('input/E_coli_v4_Build_6.experiment_feature_descriptions'))
features_long = features_long %>% select(experiment_name, feature_name, value)

# remove any NA features (found one in above dataset)
features_long = features_long %>% filter(!is.na(feature_name))

# list of features included in datasets
unique(features_long$feature_name)
```

    ##  [1] aeration                      Ampicillin                   
    ##  [3] arabinose                     cell_density                 
    ##  [5] culture_shaking               culture_temperature          
    ##  [7] culture_type                  culture_vessel               
    ##  [9] culture_vessel_volume         culture_volume               
    ## [11] experimenter                  Norfloxacin                  
    ## [13] note                          peptone                      
    ## [15] perturbation                  perturbation_gene            
    ## [17] plasmid                       resistant_to                 
    ## [19] RNA_prep_type                 RNA_stop_solution            
    ## [21] sodium_chloride               strain                       
    ## [23] structured_metadata           yeast_extract                
    ## [25] note2                         time_point                   
    ## [27] ammonium_chloride             calcium_chloride             
    ## [29] disodium_hydrogen_phosphate   glucose                      
    ## [31] growth_phase                  magnesium_chloride           
    ## [33] potassium_phosphate_dibasic   ammonium_molybdate           
    ## [35] boric_acid                    cobalt_chloride              
    ## [37] culture_pH                    cupric_sulfate               
    ## [39] ferrous_sulfate               manganese(II)_chloride       
    ## [41] MOPS                          potassium_phosphate_monobasic
    ## [43] potassium_sulfate             thiamine_HCl                 
    ## [45] tricine                       zinc_sulfate                 
    ## [47] acetate                       glycerol                     
    ## [49] l-proline                     treatment_duration           
    ## [51] Ciprofloxacin                 curator_interpretation       
    ## [53] DNase_treatment               inoculum_dilution            
    ## [55] IPTG                          succinate                    
    ## [57] Kanamycin                     RNA_type                     
    ## [59] temperature_type              antifoam                     
    ## [61] culture_O2                    potassium_nitrate            
    ## [63] HOMOPIPES                     potassium_chloride           
    ## [65] tryptone                      o-phenanthroline             
    ## [67] Spectinomycin                 dimethylformamide            
    ## [69] indole                        serine_hydroxamate           
    ## [71] l-valine                      dilution_rate                
    ## [73] WTE                           gentamicin                   
    ## [75] aluminum_potassium_disulfate  ammonium_sulfate             
    ## [77] cobalt(II)_nitrate            copper(II)_chloride          
    ## [79] disodium_EDTA                 fructose                     
    ## [81] iron(II)_ammonium_sulfate     l-cysteine                   
    ## [83] magnesium_sulfate             nickel_chloride              
    ## [85] nicotinic_acid                sodium_selenite              
    ## [87] casamino_acids                cefsulodin                   
    ## [89] mecillinam                    chloramphenicol              
    ## [91] adenine                      
    ## 91 Levels: acetate adenine aeration ... zinc_sulfate

``` {.r}
# flatten feature data
feature_data = tbl_df(dcast(features_long, experiment_name ~ feature_name,
                            value.var="value"))
dim(feature_data)
```

    ## [1] 466  92

``` {.r}
# not all of the features are guaranteed to be covered for a given sample. This
# means that there will be a number of NAs in the data:
feature_data[1:5,1:5]
```

    ## Source: local data frame [5 x 5]
    ## 
    ##    experiment_name acetate adenine        aeration
    ## 1 acid_shift_10min      NA      NA assumed_aerobic
    ## 2  acid_shift_1min      NA      NA assumed_aerobic
    ## 3  acid_shift_5min      NA      NA assumed_aerobic
    ## 4    acid_shift_t0      NA      NA assumed_aerobic
    ## 5       ast_pBAD18      NA      NA assumed_aerobic
    ## Variables not shown: aluminum_potassium_disulfate (chr)

``` {.r}
# which features which have values set for every sample?
features = colnames(feature_data)
features[complete.cases(t(feature_data))]
```

    ## [1] "experiment_name" "aeration"        "experimenter"    "strain"

``` {.r}
# distribution of strains
head(sort(table(feature_data$strain), decreasing=TRUE), 10)
```

    ## 
    ##                  MG1655             MG1655_yale                 BW25113 
    ##                     149                      75                      40 
    ##                    EMG2                     K12 MG1063 (recA56 = recA-) 
    ##                      18                      18                      16 
    ##        BW25113deltaRecA                   W3110                  MG1063 
    ##                      15                      12                      10 
    ##              BL21 (DE3) 
    ##                       8

This gives us a pretty good place to start. let's build a simplified design matrix including sample id, experiment (\~batch), and XXXX feature.

#### Limit to K12 MG1655 strain

For simplicity, let's focus for now on one strain:

-   [E. coli K12 MG1655](http://www.genome.wisc.edu/resources/strains.htm)

This is the most widely studied strain in the dataset and includes nearly half of all of the samples included in M<sup>3D</sup>.

``` {.r}
feature_data = feature_data %>% filter(strain %in% c('MG1655', 'MG1655_yale'))

# remove any features for which there are no longer any relevant samples
feature_data = feature_data[,features[colSums(!is.na(feature_data)) > 0]]
features = colnames(feature_data)

# number of samples remaining
nrow(feature_data)
```

    ## [1] 224

#### Explore features in dataset

Next, let's see which features are most commonly studied.

``` {.r}
# most common features
sort(colSums(!is.na(feature_data)), decreasing=TRUE)
```

    ##               experiment_name                      aeration 
    ##                           224                           224 
    ##                  experimenter                 RNA_prep_type 
    ##                           224                           224 
    ##                        strain               sodium_chloride 
    ##                           224                           221 
    ##             RNA_stop_solution           structured_metadata 
    ##                           220                           214 
    ##           culture_temperature                  culture_type 
    ##                           212                           188 
    ##                          note                culture_vessel 
    ##                           177                           165 
    ##                       peptone                 yeast_extract 
    ##                           164                           164 
    ##                culture_volume         culture_vessel_volume 
    ##                           161                           147 
    ##               culture_shaking                  cell_density 
    ##                           146                           127 
    ##                    time_point                  perturbation 
    ##                           118                            98 
    ##             perturbation_gene                  growth_phase 
    ##                            98                            86 
    ##                       plasmid                  resistant_to 
    ##                            84                            82 
    ##                    Ampicillin                   Norfloxacin 
    ##                            75                            71 
    ##                     arabinose                       glucose 
    ##                            69                            68 
    ##             ammonium_chloride              calcium_chloride 
    ##                            57                            57 
    ##            magnesium_chloride               DNase_treatment 
    ##                            57                            41 
    ##               ferrous_sulfate                  thiamine_HCl 
    ##                            41                            41 
    ##            ammonium_molybdate                    boric_acid 
    ##                            39                            39 
    ##               cobalt_chloride                cupric_sulfate 
    ##                            39                            39 
    ##        manganese(II)_chloride                          MOPS 
    ##                            39                            39 
    ## potassium_phosphate_monobasic             potassium_sulfate 
    ##                            39                            39 
    ##                       tricine                  zinc_sulfate 
    ##                            39                            39 
    ##             inoculum_dilution                    culture_pH 
    ##                            37                            33 
    ##                         note2   potassium_phosphate_dibasic 
    ##                            28                            20 
    ##   disodium_hydrogen_phosphate                          IPTG 
    ##                            18                            15 
    ##                      antifoam                    culture_O2 
    ##                            13                            13 
    ##                     Kanamycin              o-phenanthroline 
    ##                            13                             8 
    ##            treatment_duration                    cefsulodin 
    ##                             7                             6 
    ##                    mecillinam                 Spectinomycin 
    ##                             6                             6 
    ##            serine_hydroxamate                    gentamicin 
    ##                             4                             3 
    ##              temperature_type              ammonium_sulfate 
    ##                             3                             2 
    ##                casamino_acids                 Ciprofloxacin 
    ##                             2                             2 
    ##                      l-valine             magnesium_sulfate 
    ##                             2                             2 
    ##                      RNA_type                       acetate 
    ##                             2                             1 
    ##                       adenine                      glycerol 
    ##                             1                             1 
    ##                     l-proline                     succinate 
    ##                             1                             1

``` {.r}
# temperature
table(as.numeric(feature_data$culture_temperature))
```

    ## 
    ##  30  37  42  50 
    ##   3 197   9   3

``` {.r}
# time
# we would need to normalize the units before making any comparisons here
table(feature_data$time_point)
```

    ## 
    ##     0  1080   110    12   120   135  1440   150  1560   180    20   225 
    ##    16     1     3     1    11     1     4     2     1     5     8     1 
    ##    24   240   270    30   300   330    36   360    40   405  4320    48 
    ##     2     4     1    10     1     1     1     1     4     1     1     1 
    ##   480     5    -5    60   720 86400    90 
    ##     2     4     4    20     2     1     3

``` {.r}
# cell density
table(feature_data$cell_density)
```

    ## 
    ##  0.1  0.2  0.3 0.35  0.4    1 11.5 13.2  1.5  2.2 
    ##    1   22   81    8    8    1    2    1    2    1

``` {.r}
# growth phase
table(feature_data$growth_phase)
```

    ## 
    ##    biofilm  early_log   late_log        log stationary 
    ##          6         17         10         40         13

Growth phase exploration
------------------------

Let's start by looking at growth phase since there is a large number of samples with known growth phase information and the number of samples in each phase is decent.

``` {.r}
# discard samples with no growth phase information
feature_data = feature_data[!is.na(feature_data$growth_phase),]

# drop other unrelated features (for now...)
feature_data = feature_data %>% select(experiment_name, experimenter,
                                       growth_phase, strain)
rownames(feature_data) = 1:nrow(feature_data)
```

#### Sample design

``` {.r}
kable(feature_data)
```

|experiment\_name|experimenter|growth\_phase|strain|
|:---------------|:-----------|:------------|:-----|
|bform\_biofilm\_attachment|Ito A|biofilm|MG1655|
|bform\_biofilm\_colony|Ito A|biofilm|MG1655|
|bform\_biofilm\_maturation|Ito A|biofilm|MG1655|
|bform\_plank\_exp|Ito A|log|MG1655|
|bform\_plank\_stat|Ito A|stationary|MG1655|
|har\_S0\_noIPTG|Haddadin FT|late\_log|MG1655|
|har\_S0\_R\_noIPTG|Haddadin FT|late\_log|MG1655|
|har\_S1\_IPTG|Haddadin FT|late\_log|MG1655|
|har\_S1\_noIPTG|Haddadin FT|late\_log|MG1655|
|har\_S1\_R\_IPTG|Haddadin FT|late\_log|MG1655|
|har\_S1\_R\_noIPTG|Haddadin FT|late\_log|MG1655|
|har\_S4\_IPTG|Haddadin FT|stationary|MG1655|
|har\_S4\_noIPTG|Haddadin FT|stationary|MG1655|
|har\_S4\_R\_IPTG|Haddadin FT|stationary|MG1655|
|har\_S4\_R\_noIPTG|Haddadin FT|stationary|MG1655|
|high\_density\_heat\_shock|Haddadin FT|late\_log|MG1655|
|high\_density\_heat\_shock\_IPTG|Haddadin FT|late\_log|MG1655|
|high\_density\_serine\_hydrox\_t1|Haddadin FT|late\_log|MG1655|
|M9\_K\_appY|Covert MW|log|MG1655|
|M9\_K\_appY\_anaerobic|Covert MW|log|MG1655|
|M9\_K\_arcA|Covert MW|log|MG1655|
|M9\_K\_arcA\_anaerobic|Covert MW|log|MG1655|
|M9\_K\_arcAfnr|Covert MW|log|MG1655|
|M9\_K\_arcAfnr\_anaerobic|Covert MW|log|MG1655|
|M9\_K\_fnr|Covert MW|log|MG1655|
|M9\_K\_fnr\_anaerobic|Covert MW|log|MG1655|
|M9\_K\_oxyR|Covert MW|log|MG1655|
|M9\_K\_oxyR\_anaerobic|Covert MW|log|MG1655|
|M9\_K\_soxS|Covert MW|log|MG1655|
|M9\_K\_soxS\_anaerobic|Covert MW|log|MG1655|
|M9\_WT|Covert MW|log|MG1655|
|M9\_WT\_anaerobic|Covert MW|log|MG1655|
|MG1063\_uninduced\_t180|Kohanski MA|early\_log|MG1655|
|MG1655\_ampicillin\_t120|Kohanski MA|early\_log|MG1655|
|MG1655\_ampicillin\_t30|Kohanski MA|early\_log|MG1655|
|MG1655\_ampicillin\_t60|Kohanski MA|early\_log|MG1655|
|MG1655\_kanamycin\_t120|Kohanski MA|early\_log|MG1655|
|MG1655\_kanamycin\_t30|Kohanski MA|early\_log|MG1655|
|MG1655\_kanamycin\_t60|Kohanski MA|early\_log|MG1655|
|MG1655\_norfloxacin\_t120|Kohanski MA|early\_log|MG1655|
|MG1655\_norfloxacin\_t30|Kohanski MA|early\_log|MG1655|
|MG1655\_norfloxacin\_t60|Kohanski MA|early\_log|MG1655|
|MG1655\_spectinomycin\_t120|Kohanski MA|early\_log|MG1655|
|MG1655\_spectinomycin\_t30|Kohanski MA|early\_log|MG1655|
|MG1655\_spectinomycin\_t60|Kohanski MA|early\_log|MG1655|
|MG1655\_uninduced\_t0|Kohanski MA|early\_log|MG1655|
|MG1655\_uninduced\_t120|Kohanski MA|early\_log|MG1655|
|MG1655\_uninduced\_t30|Kohanski MA|early\_log|MG1655|
|MG1655\_uninduced\_t60|Kohanski MA|early\_log|MG1655|
|MG1655\_wt\_24hr\_biofilm|Yang X|biofilm|MG1655|
|MG1655\_wt\_R1drd19\_24hr\_biofilm|Yang X|biofilm|MG1655|
|MOPS\_K\_crp|Allen TE|log|MG1655|
|MOPS\_K\_cspA|Allen TE|log|MG1655|
|MOPS\_K\_dps|Allen TE|log|MG1655|
|MOPS\_K\_dps\_stationary|Allen TE|stationary|MG1655|
|MOPS\_K\_dps\_stationary2|Allen TE|stationary|MG1655|
|MOPS\_K\_hns|Allen TE|log|MG1655|
|MOPS\_K\_hupB|Allen TE|log|MG1655|
|rb\_wt\_biofilm|Ito A|biofilm|MG1655|
|rb\_wt\_exponential|Ito A|log|MG1655|
|rb\_wt\_stationary|Ito A|stationary|MG1655|
|rpoS\_min\_wt\_exp|Dong T|log|MG1655|
|rpoS\_min\_wt\_stat|Dong T|stationary|MG1655|
|sg\_baseline\_t0|Kohanski M|log|MG1655|
|sg\_gent\_t120|Kohanski M|log|MG1655|
|sg\_gent\_t30|Kohanski M|log|MG1655|
|sg\_gent\_t60|Kohanski M|log|MG1655|
|sg\_spect\_t120|Kohanski M|log|MG1655|
|sg\_spect\_t30|Kohanski M|log|MG1655|
|sg\_spect\_t60|Kohanski M|log|MG1655|
|sg\_untreated\_t120|Kohanski M|log|MG1655|
|sg\_untreated\_t30|Kohanski M|log|MG1655|
|sg\_untreated\_t60|Kohanski M|log|MG1655|
|WT\_MOPS\_acetate|Allen TE|log|MG1655|
|WT\_MOPS\_acidShock|Allen TE|log|MG1655|
|WT\_MOPS\_cipro|Allen TE|log|MG1655|
|WT\_MOPS\_cipro2|Allen TE|log|MG1655|
|WT\_MOPS\_glucose|Allen TE|log|MG1655|
|WT\_MOPS\_glycerol|Allen TE|log|MG1655|
|WT\_MOPS\_heatShock|Allen TE|log|MG1655|
|WT\_MOPS\_lateLog|Allen TE|late\_log|MG1655|
|WT\_MOPS\_proline|Allen TE|log|MG1655|
|WT\_MOPS\_stationary|Allen TE|stationary|MG1655|
|WT\_MOPS\_stationary2|Allen TE|stationary|MG1655|
|WT\_MOPS\_stationary3|Allen TE|stationary|MG1655|
|WT\_MOPS\_stationary4|Allen TE|stationary|MG1655|

``` {.r}
# load experiment data (includes experiment/replicate mapping)
experiment_data = tbl_df(read.delim('input/E_coli_v4_Build_6.experiment_descriptions'))
dim(experiment_data)
```

[1] 907 3

``` {.r}
# drop chips that are not related to our design
experiment_data = experiment_data %>% filter(experiment_name %in%
                                             feature_data$experiment_name)
dim(experiment_data)
```

[1] 184 3

``` {.r}
# design
design = data.frame(
    condition=feature_data$growth_phase[match(experiment_data$experiment_name,
                                              feature_data$experiment_name)],
    batch=experiment_data$experiment_name
)
```

That should be good for now. Now let's move onto the actual expression data.

### Load expression data

``` {.r}
# expression data
raw_data = tbl_df(read.delim('input/E_coli_v4_Build_6_chips907probes4297.tab.gz',
                             row.names=1))
dim(raw_data)
```

    ## [1] 4297  907

``` {.r}
raw_data = raw_data[,colnames(raw_data) %in% experiment_data$chip_name]
dim(raw_data)
```

    ## [1] 4297  184

### Helper functions

Before continuing, let's first load a few helper functions that will be useful in the downstream analysis...

``` {.r}
#
# Choose colors to use when plotting sample condition and batch
#
sample_plot_colors = function (condition, batch) {
    # Convert to factor if not already and remove any unused levels
    condition = factor(condition)
    batch = factor(batch)

    # Batch colors
    if (nlevels(batch) > 1) {
        if (nlevels(batch) <= 12) {
            rc=brewer.pal(12, "Set3")[as.integer(batch)]
        } else {
            rc=rainbow(nlevels(batch))[as.integer(batch)]
        }
    } else {
        rc = rep("green", length(batch))
    }

    # Condition colors
    if (nlevels(condition) > 1) {
        if (nlevels(condition) <= 9) {
            cc=brewer.pal(9,"Set1")[as.integer(condition)]
        } else {
            cc=tail(rainbow(nlevels(condition) +
                            nlevels(batch)),
                    nlevels(condition))[as.integer(condition)]
        }
    } else {
        cc =  rep("red",length(condition))
    }

    return(list("batch"=rc, "condition"=cc))
}

#
# Plot sample heatmap
#
plot_sample_heatmap = function (counts, condition, batch,
                                metric='dist', col='heat.colors') {
    # Compute euclidean distance or pearson correlation between samples
    if (metric == 'dist') {
        dists = dist(t(counts))
        mat = as.matrix( dists )
    } else if (metric == 'pearson') {
        mat = cor(counts)
    }

    # Select plot colors
    plot_colors = sample_plot_colors(condition, batch)

    # Heatmap plot
    hv = heatmap.2(mat, margin=c(6, 6), trace="none", key=FALSE, col=col,
                   RowSideColors=plot_colors$batch,
                   ColSideColors=plot_colors$condition)
    legend(x="topleft", legend=unique(condition),
           col=unique(plot_colors$condition), pch=15)
    #legend(x="topright", legend=unique(batch),
    #       col=unique(plot_colors$batch), pch=15)
}

#
# Plot sample PCA components
#
plot_sample_pca = function(counts, condition, batch, main="", axis1=1, axis2=2,
                           include_table=TRUE) {
    # PCA
    pca = makeSVD(counts)
    pcVar = round((pca$d^2) / sum(pca$d^2) * 100, 2)

    # X and Y axis labels
    xl = sprintf("PC%d: %.2f%% variance", axis1, pcVar[axis1])
    yl = sprintf("PC%d: %.2f%% variance", axis2, pcVar[axis2])

    # Create combined data frame
    pcaData = data.frame(SampleID=colnames(counts),
                         PC1=pca$v[,axis1], PC2=pca$v[,axis2],
                         Condition=condition, Batch=batch)

    # Plot specified principle components
    plt = ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
        geom_point(stat="identity",size=5) +
        #geom_text(aes(label=SampleID), angle=45, size=4,vjust=2) +
        xlab(xl) + ylab(yl) +
        ggtitle(sprintf("%s (PC%d vs. PC%d)", main, axis1, axis2)) +
        theme(axis.ticks=element_blank(), axis.text.x=element_text(angle=-90))
    print(plt)

    # Compute variance of each PC and how they correlate with batch and
    # condition
    if (nlevels(batch) > 1) {
        pcs = pcRes(pca$v, pca$d, condition, batch)
    } else {
        pcs = pcRes(pca$v, pca$d, condition)
    }
    rownames(pcs) = paste0("PC", rownames(pcs))

    kable(head(pcs, 30))
}
```

### Samples

``` {.r}
plot_sample_pca(raw_data,
                condition=design$condition,
                batch=design$batch)
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 86. Consider specifying shapes manually. if you must have them.
    ## Warning: Removed 166 rows containing missing values (geom_point).
    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 86. Consider specifying shapes manually. if you must have them.

![plot of chunk sample\_relationships](./README_files/figure-markdown_github/sample_relationships1.png)

||propVar|cumPropVar|cond.R2|batch.R2|
|:--|------:|---------:|------:|-------:|
|PC1|30.95|30.95|69.23|99.50|
|PC2|13.21|44.16|35.35|99.31|
|PC3|9.38|53.54|19.79|99.28|
|PC4|6.68|60.22|23.07|98.78|
|PC5|5.53|65.75|18.06|98.68|
|PC6|4.07|69.82|16.45|98.62|
|PC7|3.43|73.25|14.35|95.11|
|PC8|2.81|76.06|14.80|98.25|
|PC9|2.27|78.33|1.44|97.94|
|PC10|1.59|79.92|1.05|96.26|
|PC11|1.43|81.35|9.47|96.52|
|PC12|1.32|82.67|3.16|79.48|
|PC13|1.11|83.78|4.17|96.59|
|PC14|0.92|84.70|0.80|94.09|
|PC15|0.85|85.55|4.83|84.65|
|PC16|0.76|86.31|8.59|95.78|
|PC17|0.67|86.98|6.44|91.54|
|PC18|0.63|87.61|4.21|86.94|
|PC19|0.61|88.22|9.90|95.83|
|PC20|0.56|88.78|5.90|87.25|
|PC21|0.52|89.30|2.26|93.45|
|PC22|0.49|89.79|8.94|82.47|
|PC23|0.43|90.22|12.79|91.96|
|PC24|0.41|90.63|0.66|89.57|
|PC25|0.38|91.01|3.40|88.40|
|PC26|0.37|91.38|3.48|79.41|
|PC27|0.34|91.72|3.92|92.55|
|PC28|0.33|92.05|3.92|94.48|
|PC29|0.31|92.36|1.93|92.11|
|PC30|0.29|92.65|1.57|88.22|

``` {.r}
plot_sample_heatmap(raw_data, design$condition, design$batch)
```

![plot of chunk sample\_relationships](./README_files/figure-markdown_github/sample_relationships2.png)

### Removing batch from the data

``` {.r}
batch = design$batch

model_batch = model.matrix(~batch)
#voom_batch  = voom(normed_counts, model_batch)
#fit_batch = lmFit(voom_batch)

# Get the residuals (everything but batch effect)
#batch_residuals = residuals(fit_batch, voom_batch)

#plot_sample_pca(batch_residuals[,include], design_final, "Batch included in linear model")

#plot_sample_heatmap(batch_residuals[,include], design_final)
```

Discussion
==========

Questions
---------

1.  Dealing with batch
    -   In the context of a combined dataset such as this, does the use of `experiment_name` as batch make the most sense?
    -   Are there other better ways to divide up the dataset or to choose experiments to exclude? (e.g. remove all perturbations)
    -   Would it be possible to include more than one non-biological variable in the model, e.g. `experimenter_name`?

2.  Exploratory data analysis
    -   Any other methods for exploring the initial dataset that would be worth trying out?
        -   biplots, biological effect residuals, etc.

3.  Higher-dimension dataset
    -   In the above analysis, I started by removing all of the features, save for one batch and one biological variable of interest. Are there other ways to make use of the dataset as-is, or at least, with more of the features included?

4.  Network analysis
    -   Would this dataset be appropriate to use for constructing a co-expression or gene-regulatory network? What other information might we need?
    -   If so, what part(s) of the data should be used? and what should be excluded?

Where to next?
--------------

Some things that might be worth pursuing from here:

-   Pull in some additional annotation data:
    -   [BSgenome.Ecoli.NCBI.20080805](http://www.bioconductor.org/packages/2.13/data/annotation/html/BSgenome.Ecoli.NCBI.20080805.html), or,
    -   [org.EcK12.eg.db](http://www.bioconductor.org/packages/2.13/data/annotation/html/org.EcK12.eg.db.html)
-   Start with raw data and clean/normalize ourselves.
-   Make use of one of the larger M<sup>3D</sup> datasets, such as the one which includes intergenic regions.
-   Investigate other sources of information and data:
    -   [Reactome](http://www.reactome.org/)
    -   [EcoCyc](http://ecocyc.org/)
-   Look into more recent data sets.
-   Try integrating some RNA-Seq datasets
-   Read/discuss a paper where a E. coli network is constructed and validated using some of these resources.
-   Start trying out some of the methods we come across on our cleaned up dataset.
-   **Come up with a researh question...**

System Information
------------------

``` {.r}
sessionInfo()
```

    ## R version 3.1.0 (2014-04-10)
    ## Platform: x86_64-unknown-linux-gnu (64-bit)
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitcitations_0.6-2   RefManageR_0.8.2      bibtex_0.3-6         
    ##  [4] ggplot2_1.0.0         gplots_2.14.0         dplyr_0.2            
    ##  [7] RColorBrewer_1.0-5    reshape2_1.4          cbcbSEQ_0.9.1        
    ## [10] sva_3.11.3            genefilter_1.47.6     mgcv_1.8-0           
    ## [13] nlme_3.1-117          preprocessCore_1.27.1 corpcor_1.6.6        
    ## [16] limma_3.21.10         knitr_1.6             rmarkdown_0.2.49     
    ## [19] knitrBootstrap_1.0.0  setwidth_1.0-3        colorout_1.0-3       
    ## [22] vimcom_1.0-0         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] annotate_1.43.4      AnnotationDbi_1.27.8 assertthat_0.1      
    ##  [4] Biobase_2.25.0       BiocGenerics_0.11.3  bitops_1.0-6        
    ##  [7] caTools_1.17         colorspace_1.2-4     DBI_0.2-7           
    ## [10] digest_0.6.4         evaluate_0.5.5       formatR_0.10        
    ## [13] gdata_2.13.3         GenomeInfoDb_1.1.10  grid_3.1.0          
    ## [16] gtable_0.1.2         gtools_3.4.1         htmltools_0.2.4     
    ## [19] httr_0.3             IRanges_1.99.18      KernSmooth_2.23-12  
    ## [22] labeling_0.2         lattice_0.20-29      lubridate_1.3.3     
    ## [25] magrittr_1.0.1       markdown_0.7         MASS_7.3-33         
    ## [28] Matrix_1.1-4         memoise_0.2.1        mime_0.1.1          
    ## [31] munsell_0.4.2        parallel_3.1.0       plyr_1.8.1          
    ## [34] proto_0.3-10         Rcpp_0.11.2          RCurl_1.95-4.1      
    ## [37] RJSONIO_1.2-0.2      RSQLite_0.11.4       S4Vectors_0.1.2     
    ## [40] scales_0.2.4         splines_3.1.0        stats4_3.1.0        
    ## [43] stringr_0.6.2        survival_2.37-7      tools_3.1.0         
    ## [46] XML_3.98-1.1         xtable_1.7-3         yaml_2.1.13

``` {.r}
date()
```

    ## [1] "Thu Jul 17 09:52:45 2014"

References
----------

-   J. J. Faith, M. E. Driscoll, V. A. Fusaro, E. J. Cosgrove, B. Hayete, F. S. Juhn, S. J. Schneider, T. S. Gardner, (2007) Many Microbe Microarrays Database: Uniformly Normalized Affymetrix Compendia With Structured Experimental Metadata. <em>Nucleic Acids Research</em> <strong>36</strong> D866-D870 <a href="http://dx.doi.org/10.1093/nar/gkm815">10.1093/nar/gkm815</a>
