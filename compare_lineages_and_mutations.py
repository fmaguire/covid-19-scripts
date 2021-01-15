#!/usr/bin/env python

import pickle
import plotly.express as px
import argparse
import pandas as pd
from pathlib import Path


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")


def compare_metadata_files(previous, current):
    """
    Generate a report highlighting any new mutations or lineages in
    between the two dataframes
    """

    report = {}

    reported_mutations = [col for col in previous if col.endswith("mutation/deletion sets")]

    # summary of lineage and mutation set differences
    for category in ['pangolin_lineage'] + reported_mutations + ["division"]:

        # get changes in lineage and mutation set counts
        previous_count = previous[category].value_counts()
        previous_count.name = "Previous"
        current_count = current[category].value_counts()
        current_count.name = "Current"
        differences = pd.concat([previous_count, current_count], axis=1).fillna(0)
        differences = differences['Current'] - differences['Previous']
        differences = differences[differences != 0].sort_values(ascending=False)
        differences = differences.reset_index(name='Change in Genome Count')
        differences = differences.rename(columns={'index': category})

        # identify any new linages or mutation sets
        previous_set = set(previous[category].values)
        current_set = set(current[category].values)
        new_items = current_set - previous_set

        # add new to canada column to category df
        if len(new_items) > 0:
            differences.loc[differences[category].isin(new_items), 'Novel in Canada'] = "Yes"
            differences.loc[~differences[category].isin(new_items), 'Novel in Canada'] = "No"
        else:
            differences['Novel in Canada'] = "No"

        report[category] = differences

        # summary of changes in individual mutations and any novel ones
        if category.endswith('mutation/deletion sets'):
            individual = category.replace('sets', 'individual')
            previous[individual] = previous[category].str.split(',')
            current[individual] = current[category].str.split(',')

            previous_explode = previous.explode(individual)
            previous_mutations = previous_explode[individual].value_counts()
            previous_mutations.name = "Previous"
            current_explode = current.explode(individual)
            current_mutations = current_explode[individual].value_counts()
            current_mutations.name = "Current"

            differences_mutations = pd.concat([current_mutations,
                                               previous_mutations], axis=1).fillna(0)
            differences_mutations = differences_mutations['Current'] - differences_mutations['Previous']
            differences_mutations = differences_mutations[differences_mutations != 0].sort_values(ascending=False)
            differences_mutations = differences_mutations.reset_index(name='Change in Genome Count')
            differences_mutations = differences_mutations.rename(columns={'index': individual})

            new_mutations =  set(current_mutations.index) - set(previous_mutations.index)

            if len(new_mutations) > 0:
                differences_mutations.loc[differences_mutations[individual].isin(new_mutations), 'Novel in Canada'] = "Yes"
                differences_mutations.loc[~differences_mutations[individual].isin(new_mutations), 'Novel in Canada'] = "No"
            else:
                differences_mutations['Novel in Canada'] = "No"

            report[individual] = differences_mutations

    return report


def get_figure_height(number_categories: int) -> int:
    """
    Given a number of categories to plot gets an appropriate figure height.
    :param number_categories: The number of categories to plot.
    :return: A figure height to be used by plotly.
    """
    if number_categories < 10:
        return 400
    elif number_categories < 20:
        return 500
    elif number_categories < 30:
        return 600
    elif number_categories < 50:
        return 800
    elif number_categories < 75:
        return 1000
    elif number_categories < 100:
        return 1250
    else:
        return 5000


def apply_color(row, category):
    if row['Novel in Canada'] == 'Yes':
        row[category] = f"<span style='color:#FF0000'>{row[category]}</span>"
    return row[category]


def plot_category(report, category):
    report[category][category] = report[category].apply(lambda x: apply_color(x, category), axis=1)

    if len(report[category][category]) == 0:
        return None
    fig = px.bar(report[category], y=category, x='Change in Genome Count',
                 color='Novel in Canada', height=get_figure_height(len(report[category][category])))
    fig.update_layout(title=f"{category} monitoring",
                      yaxis=dict(title=category.replace('_', ' ').title(),
                                 tickmode='linear'))
    return fig


def summarise_report(report, report_title, output):
    if Path(output).exists():
        Path(output).unlink()

    with open(output, 'a') as f:
        f.write(f"<h1>{report_title}</h1>\n")
        for category in ['pangolin_lineage',
                         'division',
                         'S mutation/deletion sets',
                         'S mutation/deletion individual']:
            fig = plot_category(report, category)
            if fig:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
            else:
                f.write(f"<br>No change in {category.replace('_', ' ').title()}</br>")

        # too big to plot effectively so just dumping table into
        # sheet for now

        if len(report['all mutation/deletion individual']) > 0:
            f.write(report['all mutation/deletion individual'].sort_values(['Novel in Canada',
                                                                            'Change in Genome Count'],
                                                                           ascending=False).to_html(index=False))
        else:
            f.write(f"<br>No change in all mutation/deletion individual</br>")

        if len(report['all mutation/deletion sets']) > 0:
            f.write(report['all mutation/deletion sets'].sort_values(['Novel in Canada',
                                                                      'Change in Genome Count'],
                                                                           ascending=False).to_html(index=False))
        else:
            f.write(f"<br>No change in all mutation/deletion sets</br>")




if __name__ == "__main__":

    parser = argparse.ArgumentParser("Compare metadata with mutations and "
                                     "report any new linages or mutations")
    parser.add_argument("-p", "--previous",  type=check_file, required=True,
                        help="Prior days nextmeta file with mutations added")
    parser.add_argument("-c", "--current",  type=check_file, required=True,
                        help="Current days nextmeta file with mutations added")
    parser.add_argument("-o", "--output",  default="report.txt",
                        help="Output prefix")

    args = parser.parse_args()

    previous = pd.read_csv(args.previous, sep='\t')
    current = pd.read_csv(args.current, sep='\t')

    report_title = f"Canada {args.previous.parts[0]} vs {args.current.parts[0]}" \
                   f"<br>{len(current) - len(previous)} new genomes"
    report = compare_metadata_files(previous, current)
    summarise_report(report, report_title, args.output)
