# Author Louis Felix Nothias Feb 2021
import numpy as np
import pandas as pd
import altair as alt

# Visualize
def make_barchart(source):
    bars = alt.Chart(source).mark_bar().encode(
    x=alt.X('Count:Q', stack='zero', axis=alt.Axis(format='.5')),
    y=alt.Y('Annotation tool', sort='-x'),
    color=alt.Color('Annotation tool')
    )
    return bars

def make_barchart_rel(source):
    bars = alt.Chart(source).mark_bar().encode(
    x=alt.X('Count:Q', stack='zero', axis=alt.Axis(format='.0%'), title='Proportion'),
    y=alt.Y('Annotation tool', sort='-x'),
    color=alt.Color('Annotation tool')
    )
    return bars

# Matching of annotations
def make_barchart_match(source):
    bars = alt.Chart(source).mark_bar().encode(
    x=alt.X('Count:Q', stack='zero', axis=alt.Axis(format='.5')),
    y=alt.Y('Matching level', sort='-x'),
    color=alt.Color('Matching level')
    )
    return bars

def make_barchart_match_rel(source):
    bars = alt.Chart(source).mark_bar().encode(
    x=alt.X('Relative:Q', stack='zero', axis=alt.Axis(format='.0%'), title='Proportion'),
    y=alt.Y('Matching level', sort='-x'),
    color=alt.Color('Matching level')
    )
    return bars


def dist_plot(table,color_group,threshold):
    base = alt.Chart(table).encode(
        x=alt.X('SIR_MF_Zod_ZodiacScore:Q', scale=alt.Scale(domain=[threshold-0.05, 1])),
        y=alt.Y('GNPS_precursor mass:Q', scale=alt.Scale(domain=[100, 750])),
        color=color_group
    )

    points = base.mark_point()

    density_x = base.transform_density(
        density='SIR_MF_Zod_ZodiacScore',
        groupby=[color_group],
        steps=200,
        extent=[threshold-0.05, 1],
        counts=True,
        as_=['SIR_MF_Zod_ZodiacScore', 'density']
    ).mark_area(orient='vertical', opacity=0.7).encode(
        y='density:Q',
    ).properties(
        height=50
    )

    density_y = base.transform_density(
        density='GNPS_precursor mass',
        groupby=[color_group],
        steps=50,
        counts=True,
        extent=[100, 750],
        as_=['GNPS_precursor mass', 'density']
    ).mark_area(orient='horizontal', opacity=0.7).encode(
        x='density:Q',
    ).properties(
        width=50
    )

    return density_x & (points | density_y)
