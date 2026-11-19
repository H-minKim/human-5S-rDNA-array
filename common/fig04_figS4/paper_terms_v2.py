# -*- coding: utf-8 -*-
# The terminology dictionary for the publication figures, version 2.
# It corrects three things in the figure text:
# ① upper-case emphasis (HIGH/LOW/OFF/THAT/FIRST/ALL/...)
# 2) An axis named five different ways for the same thing (unit index, unit number, from the distal end, and so on).
# ③ internal tags (P10, canon, slot) in axis names
# The rules: on word boundaries only, longest first, keeping an initial capital but never capitalising a Greek letter.
# Anything inside $...$ is left alone, which would break mathtext.
import re

BETA = chr(946)       # β
MU = chr(956)         # μ
SIGMA = chr(963)      # σ
PHI = chr(966)        # φ

# -- 1) capitalised emphasis becomes ordinary text, on word boundaries ─────────────────────────────
# proper nouns (RNU2, DXZ4, CHM13, CpG, AluY, HPRC and so on) are not included
CAPS = [
    "HIGH", "LOW", "OFF", "ON", "THAT", "FIRST", "ALL", "OTHER", "IDENTICAL", "SIMULATED",
    "PROXIMAL", "DISTAL", "WIDTH", "MIDPOINT", "RATIO", "PREVIOUS", "SAME", "NOT",
    "ONE", "TWO", "BOTH", "EVERY", "NEVER", "ONLY", "MUST", "CANON", "RANK", "LAST",
    "ARRAY", "JUMP", "STEP",
    # further cases
    "PEAKS", "END", "INNER", "BIG", "BIMODAL", "OWN", "MIXED", "CLEAN", "WHOLE",
    "SAME WAY", "TRUE", "FALSE", "YES", "NO",
    "UNIT", "BIG", "SMALL", "EXACTLY", "HARD", "CALL",
    "POSITION", "LEVEL", "SHAPE", "COUNT", "WHOLE", "PART", "AGREE", "MOSTLY",
]

# -- 3) an over-long axis name becomes a noun phrase, with the method and the definition moved to the legend ────────────────
# A long axis name intrudes on the neighbouring cell so that they overlap at a width of 414.4 pt
# (133 distinct axis names are longer than 28 characters). Those below are the sentence-form ones.
LONG = {
    "Transcription coordinate tx = position inside the unit - 2121 "
    "(negative = upstream of the 5S gene, >119 = downstream spacer)":
        "Transcription coordinate (bp from the 5S gene start)",
    "Unit number from the end; inside each unit: ORF | ORF-miR+miR | LSau/spacer | (CCC/CCA)n":
        "Unit number from the array end",
    "Unit number from the end; inside each unit: AluY | spacer-1 | GA | spacer-2 | prom+gene":
        "Unit number from the array end",
    "Coupling to the unit's own state / (high-low difference / segment's dynamic range)":
        "Coupling to the unit state",
    "How far the first unit still is from the new state / (0 = already switched)":
        "Distance from the new state",
    "Unit number from the array end (each unit drawn in its reference frame)":
        "Unit number from the array end",
    "Unit number from the proximal end (each tick = 1 unit = 2,240 bp)":
        "Unit number from the proximal end",
    "Low-methylation segments / observed / expected (CpG-density null)":
        "Low-methylation segments, observed / expected",
    "Amplitude of the distance gradient / (β at d=1 minus β at d≥9)":
        "Distance-gradient amplitude",
    "Step between neighbouring units in one read (z, noise-scaled)":
        "Step between neighbouring units (z)",
    "Position inside the 2,240 bp unit (position inside the unit)":
        "Position inside the unit (bp)",
    "Distance between the two CpGs in the same molecule (bp, log)":
        "Distance between two CpGs (bp)",
    "Position of the two mixture components (unit methylation β)":
        "Fitted component mean β",
    "Tx (0 = gene start, 1-119 = gene, >119 = downstream spacer)":
        "Transcription coordinate (bp)",
    "Co-methylation φ(obs) - φ(null: same site, same haplotype)":
        "Co-methylation φ, observed − null",
    "Neighbour-CpG agreement / observed - site-preserving null":
        "Neighbour-CpG agreement, observed − null",
    "Fraction of read-units in that β bin (column-normalised)":
        "Fraction of read-units",
    "Purity = min(fraction of molecules ≤0.35, fraction ≥0.7)": "Purity",
    "Purity = min(fraction ≤0.35, fraction ≥0.7)": "Purity",
    "unit layer: bimodality separation of unit methylation β":
        "Unit layer: bimodality separation",
    "Units from the end (one box = one unit of one molecule)": "Unit number from the array end",
    "Unit index inside the read (local, not the array axis)": "Unit index inside the read",
    "Fraction of units with sub-unit methylation β 0.2-0.6":
        "Fraction of units with intermediate β",
    "Median largest single step within a molecule (|d β|)": "Median largest step (|Δβ|)",
    "Distance to nearest unit of the other state (units)": "Distance to the nearest opposite unit",
    "Within-haplotype IQR of transition midpoint (units)": "IQR of transition midpoint (units)",
    "Unit relative to the crossing (120 bins per unit)": "Unit relative to the crossing",
    "Change from distance 1 to 9+ (gradient amplitude)": "Distance-gradient amplitude",
    "Reads (start within 3 units of the end, ≥5 units)": "Reads",
    "Within-haplotype IQR of transition width (units)": "IQR of transition width (units)",
    "Position in the first unit of the new block (bp)": "Position in the first unit (bp)",
    "Sub-unit methylation β (one value per read-unit)": "Sub-unit methylation β",
    "That haplotype's mean methylation β, distal 1-5": "Mean β, distal units 1-5",
    "Probability-mean methylation β (ML/255 per CpG)": "Probability-mean β",
    "CpG layer: positional fixation of the low state": "CpG layer: positional fixation",
    "Share of a molecule's rise from big steps (%)": "Share of rise from large steps (%)",
    "Score (0 = low interior, 1 = high interior)": "Normalised state",
    "Position inside the unit (bp) . canon frame": "Position inside the unit (bp)",
    "Score in the first unit of a new low block": "Score in the first unit",
    "Mean unit methylation β of that chromosome": "Mean unit methylation β",
    "Unit methylation β of that single molecule": "Unit methylation β",
    "Standard deviation of segment mean $\\beta$": "SD of segment mean $\\beta$",
    "Fraction of off units (%) . complete units": "Fraction of off units (%)",
    "Fraction of off units (%) cut 0.4611": "Fraction of off units (%)",
    "Window start φ (position inside the unit)": "Window start (bp)",
    "Mean methylation β of a dispersed 5S site": "Mean methylation β of a site",
    "Observed / expected (vs CpG-density null)": "Observed / expected",
    "Mean rise per read in the transition zone": "Mean rise per read",
    "Score in the first unit of the new block": "Score in the first unit",
    "Unit number from distal end (1 = distal)": "Unit number from the distal end",
    "Mean methylation β of that compartment": "Mean compartment β",
    "Unit layer: Ashman D (>2 = two peaks)": "Unit layer: Ashman D",
    "Fraction of read-units / with β ≥ 0.5": "Fraction of read-units with β ≥ 0.5",
    "Intercept = that compartment's floor": "Compartment floor (intercept)",
    "Mean methylation β of that haplotype": "Mean haplotype β",
    "Unit methylation β of that molecule": "Unit methylation β",
    "% of the unit's common bins covered": "Common bins covered (%)",
    "Share of one molecule's total rise": "Share of the total rise",
    "Share of a molecule's total rise": "Share of the total rise",
    "Methylation of one molecule's core": "Core methylation β",
    "Methylation of one molecule's unit": "Unit methylation β",
    "Share of units with exactly 0 (%)": "Units with exactly zero (%)",
    "Mean CpG methylation (hard call)": "Mean CpG methylation β",
    "Change in β when that side flips": "Change in β on flipping",
    "Rise from small steps (per read)": "Rise from small steps",
    "Rise from big jumps (per read)": "Rise from large steps",
    "Where the 'low peak' lands (β)": "Position of the low peak (β)",
    "Coupling = observed φ - null φ": "Coupling (observed − null φ)",
    "Score of the block's first unit": "Score in the first unit",
    "Median coverage of that element": "Median coverage",
    "% of that direction's switches": "Switches in that direction (%)",
    "첫 유닛에서 아직 안 넘어간 정도 / "
    "(0 = 이미 다 넘어갔다)": "새 상태까지 남은 거리",
}
def _norm(x):
    """Normalise by removing line breaks, runs of spaces and ' / '.

    A line break written as ' / ' in a dictionary key
    does not match the real label (which has a line break).
    """
    return " ".join(x.replace(" / ", " ").split())


_LONG = sorted(((_norm(k), v) for k, v in LONG.items()), key=lambda kv: -len(kv[0]))

# -- 2) different words for the same thing become one, ignoring case ────────────────
PHRASE = [
    # the axis counting positions: five names for one thing
    (r"unit index from (?:the )?(?:dist|distal) (?:end|junction)",
     "unit number from the distal end"),
    (r"unit index from (?:the )?(?:prox|proximal) (?:end|junction)",
     "unit number from the proximal end"),
    (r"unit number counted from the (distal|proximal) end", r"unit number from the \1 end"),
    (r"unit index from that end", "unit number from the array end"),
    (r"unit number from each end", "unit number from the array end"),
    # the name of the ruler, with beta written as a symbol
    (r"\bunit beta\b", "unit methylation " + BETA),
    (r"\bread-unit beta\b", "read-unit methylation " + BETA),
    (r"\bmean beta\b", "mean methylation " + BETA),
    (r"\bmethylation beta\b", "methylation " + BETA),
    (r"\bgene-body beta\b", "gene-body methylation " + BETA),
    (r"\bsub-unit beta\b", "sub-unit methylation " + BETA),
    (r"\braw methylation beta\b", "raw methylation " + BETA),
    (r"\bbeta\b", BETA),
    # internal tags are removed from the axis
    (r"\s*[-–—.]\s*P10,? \d+ individuals", ""),
    (r"\(P10,? \d+ ind\.?\)", ""),
    (r"\bunit_pos\b", "position inside the unit"),
    (r"\bD4Z4_4q\b", "D4Z4 4q"),
    (r"\bD4Z4_10q\b", "D4Z4 10q"),
    (r"\bSST1_a\b", "SST1"),
    # a tick name long enough to run into its neighbour (in Claim_G c, J3 and J1)
    (r"\bAluY element\b", "AluY"),
    (r"\bDXZ4_F\b", "DXZ4, female"),
    (r"\bDXZ4_M\b", "DXZ4, male"),
    (r"\bDUX4_ORF\b", "DUX4 ORF"),
    (r"\bSAT-VAR\b", "satellite variant"),
    # abbreviations
    (r"\bprox\b", "proximal"),
    (r"\bdist\b", "distal"),
    (r"\bH2L\b", "high to low"),
    (r"\bL2H\b", "low to high"),
    (r"\bhigh-methyl\b", "high-methylation"),
    (r"\blow-methyl\b", "low-methylation"),
    (r"\bind\.", "individuals"),
    # a mathematical symbol is written as a symbol
    (r"\bmu\b", MU),
    (r"\bsigma\b", SIGMA),
    (r"\bphi\b", PHI),
    (r"\bfrac\b", "fraction"),
    # an inequality is written in the journal's notation
    (r">=", chr(8805)),
    (r"<=", chr(8804)),
]

_CAPS_RE = [(re.compile(r"\b" + w + r"\b"), w.lower())
            for w in sorted(CAPS, key=len, reverse=True)]
_PH_RE = [(re.compile(a, re.I), b) for a, b in PHRASE]


def _sub_outside_math(t):
    """The inside of $...$ math is left alone: mathtext would break."""
    parts, i = [], 0
    for m in re.finditer(r"\$[^$]*\$", t):
        parts.append((t[i:m.start()], True))
        parts.append((m.group(0), False))
        i = m.end()
    parts.append((t[i:], True))
    out = []
    for seg, do in parts:
        if do and seg:
            for rx, rep in _CAPS_RE:
                seg = rx.sub(rep, seg)
            for rx, rep in _PH_RE:
                seg = rx.sub(rep, seg)
        out.append(seg)
    return "".join(out)


def apply_terms2(t):
    if not t:
        return t
    head_upper = t[0].isalpha() and t[0].isupper()
    t = _sub_outside_math(t)
    # a long sentence-form axis name becomes a noun phrase, with the definition and the method moved to the legend, then matched by character
    _k = _norm(t)
    for a, b in _LONG:
        if _k == a:
            t = b
            break
        if a in _k:
            t = _k.replace(a, b)
            break
    t = re.sub(r"[ ]{2,}", " ", t).strip()
    # a Greek letter is never capitalised (beta written as the symbol would become a capital Beta)
    if head_upper and t and "a" <= t[0] <= "z":
        t = t[0].upper() + t[1:]
    return t
