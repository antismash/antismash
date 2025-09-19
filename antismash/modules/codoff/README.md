# codoff Module for antiSMASH

The codoff module provides codon usage analysis for biosynthetic gene clusters (BGCs) within antiSMASH. It assesses the similarity of codon usage patterns between a BGC and its genomic context, which can provide insights into horizontal gene transfer events and translation optimization.

## Overview

This module implements the core functionality from the [codoff repository](https://github.com/Kalan-Lab/codoff):

- **Codon frequency analysis**: Calculates codon usage frequencies for both BGC regions and genome-wide background
- **Statistical comparison**: Uses cosine distance and Spearman correlation to compare usage patterns
- **Empirical P-values**: Monte Carlo simulations to assess statistical significance
- **Multi-scaffold support**: Handles genomes with multiple scaffolds by analyzing them jointly
- **HTML visualization**: Generates histograms and detailed reports in the antiSMASH output

> [!NOTE] 
> The antiSMASH codoff module may find fewer total CDS features compared to the standalone codoff program because only records with BGC regions are considered. This should usually be valid, but be cautious of the results if many scaffolds lack BGC regions as the total codon usage of the genome might not be properly inferred. 

## Usage

### Command Line Options

The module is **optional** and must be explicitly enabled:

```bash
antismash --codoff-enabled [options] input_file
```

**Available options:**
- `--codoff-enabled`: Enable codon usage analysis
- `--codoff-simulations N`: Number of Monte Carlo simulations (default: 10,000)
- `--codoff-min-genome-size N`: Minimum genome size in bp (default: 500,000)
- `--codoff-max-genome-size N`: Maximum genome size in bp (default: 50,000,000)

### Example

```bash
antismash --codoff-enabled --codoff-simulations 5000 input.gbk
```

## Genome Size Requirements

The module has built-in genome size validation:
- **Minimum genome size**: Default is 500,000 bp
- **Maximum genome size**: Default is 50,000,000 bp
- **Multi-scaffold**: All scaffolds are considered jointly for size calculation

Genomes outside these limits will be skipped with appropriate warning messages.

## Output

### HTML Report

The module adds a new "codoff Analysis" tab to each region in the antiSMASH HTML report, containing:

1. **Important Notices**: Usage guidelines and multiple testing correction recommendations
2. **Statistical Results**: Empirical P-value, cosine distance, Spearman correlation
3. **Interpretation**: Clear explanation of significance (p < 0.05 threshold)
4. **Histogram**: Distribution of simulated cosine distances with observed value marked
5. **Codon Frequency Table**: Detailed comparison of standardized frequencies

### Statistical Significance

- **p < 0.05**: Significant difference in codon usage, potentially indicating horizontal gene transfer or atypical expression/translation frequency
- **p ≥ 0.05**: No significant difference, suggesting similar codon usage to genomic context

### Files Generated

- **HTML sections**: Integrated into the main antiSMASH report
- **Histogram plots**: SVG files in `codoff_plots/` subdirectory
- **JSON results**: Stored in the main antiSMASH results file

## Technical Details

### Algorithm

The module uses a **deferred analysis approach** to ensure statistical consistency across multi-scaffold genomes:

1. **Data Collection Phase**: During `run_on_record()` for each scaffold:
   - Gathers CDS features and BGC regions from the current scaffold
   - Accumulates data in global state (`_codoff_global_state`)
   - Stores BGC region metadata for later analysis

2. **Analysis Trigger**: During antiSMASH's `write_outputs()` phase:
   - **Codon Counting**: Extracts nucleotide sequences and counts codons in 3-base windows
   - **Frequency Calculation**: Computes standardized frequencies (proportions) for comparison
   - **Statistical Analysis**: Calculates cosine distance and Spearman correlation
   - **Monte Carlo Simulation**: Uses serial processing to randomly sample CDS features and generate null distribution (matches standalone codoff)
   - **P-value Calculation**: Empirical significance based on simulation results

### Multi-scaffold Handling

The module uses a **global state system** to:
- Maintain CDS feature mappings to their source records for accurate sequence extraction
- Ensure consistent statistical comparisons by analyzing all BGCs against the complete genomic background

### Complex CDS Feature Support

The module properly handles complex CDS features including:
- **CDS with introns**: Automatically extracts spliced (mature) sequences
- **Multiple exons**: Joins exonic sequences correctly
- **Reverse complement genes**: Handles minus-strand features automatically
- **Ambiguous nucleotides**: Filters out codons containing N, R, Y, etc.

This is achieved using BioPython's `location.extract()` method which handles all splicing operations.

### Performance Considerations

- **Deferred analysis**: Computational work is delayed until antiSMASH's post-processing phase, avoiding redundant calculations
- **Serial processing**: Monte Carlo simulations run sequentially to match standalone codoff behavior and ensure identical results
- **Simulation count**: Higher values provide more accurate P-values but increase runtime
- **CDS filtering**: Only processes CDS features with length divisible by 3 (matches standalone codoff)  
- **Memory usage**: Efficient storage of CDS-to-record mappings in global state
- **Single execution**: Analysis runs only once per antiSMASH run via the `write_outputs()` mechanism

## Module Structure

The codoff module follows antiSMASH's standard module patterns:

- **Clean function docstrings**: Brief, consistent documentation style
- **Minimal logging**: Essential warnings and errors only, no verbose debug output
- **Consistent naming**: All module references use lowercase "codoff"
- **Standardized test format**: Test docstrings and structure match antiSMASH conventions

## Dependencies

**Required Python packages:**
- `numpy`: Numerical operations and array handling
- `scipy`: Statistical calculations (correlation, distance metrics)
- `matplotlib`: Plot generation

**antiSMASH integration:**
- `antismash.common.secmet`: Core data structures
- `antismash.common.html_renderer`: HTML output generation
- `antismash.config`: Configuration management

## Testing

Run the module tests with:

```bash
python -m pytest antismash/modules/codoff/test/ -v
```

The test suite covers:
- Core analysis functions
- Statistical calculations (serial processing)
- Results serialization
- Deferred analysis workflow via write_outputs()
- Mock data handling
- Algorithm accuracy validation
- Integration tests with real data files

## Limitations

- **CDS length filtering**: Only requires CDS length divisible by 3 (matches standalone codoff)
- **Codon extraction**: Assumes standard 3-base codon structure, coding genes not divisible by 3 are excluded
- **Memory usage**: Large genomes with many CDS features may require significant memory