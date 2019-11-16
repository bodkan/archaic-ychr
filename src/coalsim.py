#!/usr/bin/env python3

import argparse

import pandas as pd
import msprime as msp


GEN_TIME = 25


def years_to_gen(y):
    return int(y  / GEN_TIME)


def get_id(sample, pop_params):
    for i, ind in enumerate(all_inds(pop_params)):
        if ind == sample:
            return i

def get_tmrca(ts, ind1, ind2, pop_params):
    return ts.first().tmrca(get_id(ind1, pop_params), get_id(ind2, pop_params)) * GEN_TIME


def define_samples(pop_params):
    """Generate list of sample definitions for msprime."""
    sample_names = []
    for i, pop in enumerate(pop_params):
        times = [years_to_gen(t) for t in pop_params[pop]["t_sample"]]
        sample_names.extend([msp.Sample(population=i, time=t) for t in times])
    return sample_names


def all_inds(pop_params):
    """Generate list of all simulated sample names."""
    sample_names = []
    for p in pop_params:
        n_pop = len(pop_params[p]["t_sample"])
        sample_names.extend([f"{p}{i}" for i in range(n_pop)])
    return sample_names


def Ne_scale(Ne, chrom, R):
    """Convert Ne to Y/mtDNA scale and scale by female-to-male ratio."""
    if chrom == "mt":
        Ne_scaled = Ne / (1 + 1 / R)
    elif chrom == "Y":
        Ne_scaled =  Ne / (1 + R)
    return int(1/2 * Ne_scaled) # scaling due to msprime's default diploid model


parser = argparse.ArgumentParser()

parser.add_argument("--Ne-nonafr", type=int, default=5000)
parser.add_argument("--Ne-afr", type=int, default=10000)
parser.add_argument("--Ne-arch", type=int, default=1000)

parser.add_argument("--chrom", choices = ["mt", "Y"], required=True)
parser.add_argument("--female-to-male", type=float, default=1.0)

parser.add_argument("--seq-len", help="Sequence length", type=int, default=7_000_000)
parser.add_argument("--mut-rate", help="Mutation rate [per bp per year]",
    type=float, default=7e-10)

parser.add_argument("--ui-age", nargs="+", type=int,
    help="Age of the Ust-Ishim individual [years BP]")
parser.add_argument("--arch-ages", nargs="+", type=int,
    help="Ages of archaic samples [years BP]", required=True)
parser.add_argument("--neur", type=int, help="Number of European chromosomes to simulate",
    required=True)
parser.add_argument("--nafr", type=int, default=1,
    help="Number of African chromosomes to simulate")

parser.add_argument("--split-chimp", type=int, default=6_000_000)
parser.add_argument("--split-arch", type=int, required=True)
parser.add_argument("--split-afr", type=int, required=True)

parser.add_argument("--output", metavar="FILE", help="Output filename")
parser.add_argument("--format", choices=["snp", "fa"], required=True)

parser.add_argument("--debug", action="store_true", help="Print debugging info",
    default=False)

args = parser.parse_args()
#args = parser.parse_args("--debug --arch-ages 130000 --split-arch 650_000 --split-afr 250_000 "
# "--ui-age 45000 --neur 50 --nafr 1 --output out --format snp".split())

pop_params = {
    "chimp": {"id": 0, "Ne": Ne_scale(args.Ne_afr, args.chrom, args.female_to_male),
        "t_sample": 1 * [0], "t_split" : args.split_chimp},
    "arch": {"id": 1, "Ne": Ne_scale(args.Ne_arch, args.chrom, args.female_to_male),
        "t_sample": args.arch_ages, "t_split": args.split_arch},
    "afr": {"id": 2, "Ne": Ne_scale(args.Ne_afr, args.chrom, args.female_to_male),
        "t_sample": args.nafr * [0]},
    "eur": {"id": 3, "Ne": Ne_scale(args.Ne_nonafr, args.chrom, args.female_to_male),
        "t_sample": args.ui_age + args.neur * [0], "t_split": args.split_afr},
}

id_chimp, id_arch, id_afr, id_eur  = [pop_params[p]["id"] for p in pop_params]

t_chimp, t_arch, t_eur  = [
    years_to_gen(pop_params[p]["t_split"])
    for p in ["chimp", "arch", "eur"]
]

demography = [

    # out of Africa migration
    msp.MassMigration(t_eur, id_eur, id_afr, 1.0),

    # split of the archaic lineage from modern humans
    msp.MassMigration(t_arch, id_arch, id_afr, 1.0),

    # split of chimps from the human lineage
    msp.MassMigration(t_chimp, id_chimp, id_afr, 1.0),
]

demography.sort(key=lambda event: event.time)

pop_config = [msp.PopulationConfiguration(initial_size=pop_params[p]["Ne"]) for p in pop_params]
samples = define_samples(pop_params)

if args.debug:
    print("Population identifiers:")
    for p, i in pop_params.items():
        print(f"    {p}: {i['id']}")
    msp.DemographyDebugger(
        population_configurations=pop_config,
        demographic_events=demography
    ).print_history()

# enforce the phylogeny (chimp; (archaic; (african; (eur1, eur2)))) by
# restarting the simulation if necessary
while True:
    ts = msp.simulate(
        length=args.seq_len,
        mutation_rate=args.mut_rate * GEN_TIME,
        recombination_rate=0,
        samples=samples,
        population_configurations=pop_config,
        demographic_events=demography
    )
    if get_tmrca(ts, "chimp0", "arch0", pop_params) > get_tmrca(ts, "arch0", "afr0", pop_params) and \
       all(get_tmrca(ts, "arch0", f"eur{i}", pop_params) > get_tmrca(ts, "afr0", f"eur{j}", pop_params) for i in range(args.neur + 1) for j in range(args.neur + 1)):
         break
    else:
         print(f"replicate {args.output} failed - restarting...")

snps = pd.DataFrame(
    ts.genotype_matrix(),
    columns=all_inds(pop_params),
    index=(v.position for v in ts.variants())
).rename(columns={"chimp0" : "chimp", "eur0": f"ustishim"})

snps.to_csv(args.output + ".snp.gz", sep="\t", index_label="pos")
with open(args.output + ".fa", "w") as fa:
    for s in snps.columns:
        print(f">{s}", file=fa)
        print("".join("A" if i == 0 else "T" for i in snps[s]), file=fa)
with open(args.output + ".txt", "w") as tmrca:
    print("tmrca_afr\ttmrca_arch", file=tmrca)
    print(f"{get_tmrca(ts, 'afr0', 'eur1', pop_params)}\t{get_tmrca(ts, 'arch0', 'eur1', pop_params)}", file=tmrca)
