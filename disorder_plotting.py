import pandas
import numpy
from matplotlib import pyplot
import scipy.stats
import matplotlib.patches


# Class to allow creation of different plots and potentially statistical analysis from provided list of proteins and prediction data
class DisorderAnalysis:
    # For initialisation:
    #   At least one protein list (csv) stored in Data/Protein_lists
    #   A dataframe containing iupred prediction of all proteins in lists 
    #   A disorder_threshold, defeault is set to 0.5
    #   A region_size: default is 10
    def __init__(
        self, protein_lists, prediction_data, disorder_threshold=0.5
    ):
        self.prediction_df = prediction_data

        self.protein_dict = {}
        self.disordered_residues_dict = {}
        self.sequence_lengths_dict = {}
        self.disorder_ratio_dict = {}
        self.disorder_regions_dict = {}
        for key in protein_lists:
            # Reads in specified protein lists from Data/Protein_lists/
            file_path = f"Data/Protein_lists/{key}.csv"
            self.protein_dict[key] = pandas.read_csv(
                file_path, index_col=0, names=["Protein"]
            )

            # Calculates which sequence elements meet the disorder threshold
            data = self.prediction_df[
                self.prediction_df.columns.intersection(
                    self.protein_dict[key]["Protein"]
                )
            ]
            disordered_residues = data.ge(disorder_threshold, axis="columns")

            # Calculates the ratio of disordered residues compared to the total sequence for each protein
            disorder_sum = disordered_residues.sum()
            sequence_lengths = data.ge(0, axis="columns").sum()
            disorder_ratio = disorder_sum / sequence_lengths

            # Saves values as class variables
            self.disordered_residues_dict[key] = disordered_residues
            self.sequence_lengths_dict[key] = sequence_lengths
            self.disorder_ratio_dict[key] = disorder_ratio
            # self.disorder_regions_dict[key] = pandas.Series(disordered_regions_series)

    def region_counter(self, protein_keys, region_sizes):

        idr_stats_dict = {}
        for key in protein_keys:
            # Retrieve disordered residues data for the current key
            disordered_residues = self.disordered_residues_dict[key]
            protein_stats_dict = {}

            # Iterate over each protein in the disordered residues data
            for protein in disordered_residues.columns:
                # Identify transitions (changes) in disordered residues
                sequence_array = numpy.array(disordered_residues[protein])
                transitions = numpy.diff(
                    numpy.concatenate(([False], sequence_array, [False]))
                )

                # Identify starts of disordered regions
                region_starts = numpy.where(transitions == True)[0]
                # Calculate lengths of disordered regions
                disordered_regions = numpy.array(
                    region_starts[1::2] - region_starts[0::2]
                )

                # Calculate protein statistics based on disordered regions
                if numpy.size(disordered_regions) == 0:
                    protein_stats_dict[protein] = [0, 0, 0]
                else:
                    protein_stats_dict[protein] = [
                        disordered_regions.sum(),
                        numpy.size(disordered_regions),
                        disordered_regions.max(),
                    ]

            # Create DataFrame from protein statistics dictionary
            idr_stats_df = pandas.DataFrame.from_dict(
                protein_stats_dict,
                orient="index",
                columns=["Disordered_residues", "Region_count", "Largest_region"],
            )
            idr_stats_dict[key] = idr_stats_df

        idr_percentages_dict = {}
        # Iterate over each region size
        for size in region_sizes:
            # Calculate percentage of proteins with largest disordered region >= current region size
            protein_count = numpy.array(
                [len(idr_stats_dict[key]["Largest_region"]) for key in protein_keys]
            )
            sum_disordered = numpy.array(
                [
                    idr_stats_dict[key]["Largest_region"].ge(size).sum()
                    for key in protein_keys
                ]
            )
            protein_ratio = sum_disordered / protein_count * 100

            idr_percentages_dict[size] = protein_ratio

        # Create DataFrame from the percentage dictionary, with protein keys as index
        idr_percentages_df = pandas.DataFrame(idr_percentages_dict, index=protein_keys)

        return idr_percentages_df

    def bootstrap_error(self, protein_key, region_sizes, resamples=100):

        disordered_residues = self.disordered_residues_dict

        bootstraps = []
        for region_size in region_sizes:

            def largest_region(sequence_array):
                # Identify transitions (changes) in disordered residues
                transitions = numpy.diff(
                    numpy.concatenate(([False], sequence_array, [False]))
                )

                # Identify starts of disordered regions
                region_starts = numpy.where(transitions == True)[0]
                # Calculate lengths of disordered regions
                disordered_regions = numpy.array(
                    region_starts[1::2] - region_starts[0::2]
                )

                # if numpy.size(disordered_regions) == 0:
                #     region_max = 0
                # else:
                #     region_max = disordered_regions.max()

                return (
                    disordered_regions.max()
                    if numpy.size(disordered_regions) > 0
                    else 0
                )

            def statistic(sample):
                protein_count = len(sample)
                sum_disordered = numpy.sum(sample >= region_size)

                return (sum_disordered / protein_count) * 100

            disordered_residues_array = numpy.array(disordered_residues[protein_key])
            largest_regions = numpy.apply_along_axis(
                largest_region, axis=0, arr=disordered_residues_array
            )

            bootstraps.append(
                1.96
                * scipy.stats.bootstrap(
                    (largest_regions,), statistic, n_resamples=resamples
                ).standard_error
            )

        return bootstraps
    
    def disorder_ratio_proportion(self, protein_keys):
        # Define step size for intervals
        step_size = 0.1
        distribution_dict = {}

        for key in protein_keys:
            # Retrieve disorder ratio data for the current protein key
            data = numpy.array(self.disorder_ratio_dict[key])
            protein_count = len(data)

            # Calculate distribution points based on intervals
            distribution_points = numpy.array(
                [
                    numpy.count_nonzero(
                        (data >= i * step_size) & (data <= (i + 1) * step_size)
                    )
                    / protein_count
                    for i in range(10)
                ]
            )

            # Store distribution points in the dictionary for the current protein key ( starting from zero)
            distribution_dict[key] = numpy.concatenate(([0], distribution_points))

        return pandas.DataFrame(distribution_dict)

    def bar_plot_regions(
        self, protein_keys, region_sizes, labels, resamples=100, fig_size=None
    ):

        if fig_size == None:
            fig_size = (9 + len(region_sizes), 10)

        # Get a Dark2 colormap and create a list of colors based on protein keys
        cmap = pyplot.get_cmap("Dark2")
        colors = [cmap(i) for i in range(len(protein_keys))]

        # Call region_counter to get percentage of proteins for each region size and transpose the DataFrame
        idr_percentages_df = self.region_counter(protein_keys, region_sizes).transpose()

        idr_bootstrap_error = pandas.DataFrame.from_dict(
            {
                key: (self.bootstrap_error(key, region_sizes, resamples))
                for key in protein_keys
            }
        )
        idr_bootstrap_error.index = region_sizes

        # Plot DataFrame as a bar chart
        ax = idr_percentages_df.plot(
            kind="bar",
            yerr=idr_bootstrap_error,
            color=colors,
            figsize=fig_size,
            width=0.6,
            edgecolor="#122423",
            zorder=3,
            capsize=4,
            rot=0,
        )

        ax.grid(color="grey", axis="y", zorder=0)
        ax.set_ylabel("Percent of total protein set", fontsize = 15)
        ax.set_xlabel("Longest disordered region >=", fontsize = 15)
        ax.tick_params(labelsize=14)
        ax.legend(labels, fontsize = 15)
        
        return ax

    def plot_histogram(self, protein_keys, labels):

        # Calculates percentege portions of all disorder ratios
        distribution_df = self.disorder_ratio_proportion(protein_keys)

        step_size = 0.1

        # figure size
        fig, ax = pyplot.subplots(figsize=(9.2, 5))

        # Add gridlines along the x-axis and set x-axis properties
        ax.grid(color="grey", axis="x", zorder=0)
        ax.set_xticks(numpy.linspace(0, 1, 11))
        ax.set_xlim(0, 1)
        ax.set_xlabel("Portion of total protein count")

        # Get a colormap and define category colors and names based on intervals
        category_colors = reversed(
            pyplot.get_cmap("summer")(numpy.linspace(0.15, 0.85, 5))
        )
        category_names = [
            f"{100 * step_size*2 * i} - {100 * step_size*2 *(i+1)}% " for i in range(5)
        ]

        start = 0
        split = 0

        # Iterate over each category (interval) and plot horizontal bars
        for i, (color, category) in enumerate(zip(category_colors, category_names)):
            widths = distribution_df.loc[(i * 2) + 1 : (i + 1) * 2, :].sum()
            bars = ax.barh(
                labels,
                widths,
                left=start,
                height=0.5,
                label=category,
                color=color,
                zorder=2,
                edgecolor="black",
            )
            start += widths

            # Adds lines at each other 10% interval
            split += distribution_df.loc[(i * 2) : (i * 2) + 1, :].sum()

            for bar, x in zip(bars, split):
                height = bar.get_height()
                y_bottom = bar.get_y()  # Bottom y-coordinate of the bar
                y_top = y_bottom + height  # Top y-coordinate of the bar

                ymin_normalized = (y_bottom - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
                ymax_normalized = (y_top - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
                
                # Add vertical dashed line at the center of the bar
                ax.axvline(x, ymin=ymin_normalized, ymax=ymax_normalized, color="black", linestyle="dashed", zorder=3)

        ax.legend(
            title="Disorder Ratio Intervals",
            ncol=5,
            bbox_to_anchor=(0, 1),
            loc="lower left",
            fontsize="small",
        )

        return fig, ax

    def mannwhitneyu_plot(self, protein_keys, labels=None):
        # Extract disorder data and calculate pairwise Mann-Whitney U test p-values
        disorder_data = self.disorder_ratio_dict
        p_values = []

        for key1 in protein_keys:
            for key2 in protein_keys:
                _, p_value = scipy.stats.mannwhitneyu(
                    disorder_data[key1], disorder_data[key2]
                )
                p_values.append((key1, key2, p_value))

        # Create DataFrame with calculated p-values and corrected p-values
        p_values_df = pandas.DataFrame(p_values, columns=["key1", "key2", "p_value"])
        p_values_df["p_value_corrected"] = scipy.stats.false_discovery_control(
            p_values_df["p_value"], method="by"
        )

        # Pivot table for corrected p-values
        matrix = pandas.pivot_table(
            p_values_df, values="p_value_corrected", index="key1", columns="key2"
        )
        protein_keys_rev = matrix.columns.copy()
        matrix.rename(
            columns={key: label for key, label in zip(protein_keys, labels)},
            inplace=True,
        )

        # Set up plot
        pyplot.figure(figsize=(8, 6))
        cell_size = 1.0
        colors = pyplot.get_cmap("summer")(
            numpy.linspace(0.1, 0.9, 2)
        )  # Use a coolwarm colormap for colors

        # Loop through protein keys for pairwise comparisons
        for i, key1 in enumerate(protein_keys_rev):
            for j, key2 in enumerate(protein_keys_rev):
                p_val = matrix.iloc[i, j]

                if p_val < 0.05:  # Consider only significant comparisons (alpha = 0.05)
                    median_class1 = disorder_data[key1].median()
                    median_class2 = disorder_data[key2].median()

                    if median_class1 > median_class2:
                        color = colors[0]  # Use first color for class1 > class2
                    elif median_class1 < median_class2:
                        color = colors[1]  # Use second color for class1 < class2
                    else:
                        continue  # Skip if medians are equal

                    # Plot rectangles for significant comparisons
                    rect = pyplot.Rectangle(
                        (j - 0.5, i - 0.5),
                        cell_size,
                        cell_size,
                        facecolor=color,
                        edgecolor="black",
                        zorder=0,
                    )
                    pyplot.gca().add_patch(rect)

        # Customize plot settings and labels
        pyplot.xticks(
            ticks=numpy.arange(len(matrix.columns)), labels=matrix.columns, rotation=45
        )
        pyplot.yticks(
            ticks=numpy.arange(len(matrix.columns)), labels=matrix.columns, rotation=0
        )

        # Add a diagonal line
        pyplot.plot(
            [-0.5, len(matrix.columns) - 0.5],
            [-0.5, len(matrix.index) - 0.5],
            color="black",
            linestyle="--",
            linewidth=1.0,
        )

        # Add grid lines between cells
        for i in range(len(matrix.columns)):
            pyplot.axvline(
                x=i - 0.5, color="black", linestyle="-", linewidth=0.5, zorder=3
            )
            pyplot.axhline(
                y=i - 0.5, color="black", linestyle="-", linewidth=0.5, zorder=3
            )

        pyplot.xlim(-0.5, len(matrix.columns) - 0.5)
        pyplot.ylim(
            len(matrix.index) - 0.5, -0.5
        )  # Reverse y-axis to match matrix orientation
        pyplot.tight_layout()
        
        pop_a = matplotlib.patches.Patch(color= colors[0], label='Row median greater') 
        pop_b = matplotlib.patches.Patch(color=colors[1] , label='Row median lesser') 
        
        pyplot.legend(
            title = 'Pairwise Comparison of Disorder Content Medians',
            handles = [pop_a, pop_b],          
            ncol=2,
            bbox_to_anchor=(0.25, 1),
            loc="lower left",
            fontsize="medium",)

        return pyplot.show()

    # Plots boxplot
    def plot_box(self, protein_keys, labels, fig_size = (9, 6)):
        
        data = [self.disorder_ratio_dict[key] for key in protein_keys]
        
        pyplot.figure(figsize = fig_size)
        pyplot.boxplot(data, labels = labels, showmeans=True, meanline=True)
        pyplot.title('Disorder content')
        
        return pyplot.show() 


    # Collects some interesting parameters for each protein-key in a dataframe
    def tabulated_values(self, protein_keys):
        tab_dict = {}
        for key in protein_keys:
            tab_dict[key] = numpy.around(
                numpy.array(
                    [
                        len(self.sequence_lengths_dict[key]),
                        self.sequence_lengths_dict[key].max(),
                        self.sequence_lengths_dict[key].mean(),
                        self.sequence_lengths_dict[key].median(),
                        self.disorder_ratio_dict[key].mean(),
                        self.disorder_ratio_dict[key].median(),
                    ]
                ),
                3,
            )
        tabulated_df = pandas.DataFrame(tab_dict)
        tabulated_df.index = [
            "Protein count",
            "Max sequence length",
            "Average sequence length",
            "Median sequence length",
            "Average disorder ratio",
            "median disorder ratio",
        ]
        return tabulated_df
