# Plot Standardization Action Plan

**Date:** 2025-11-05
**Status:** Ready for Implementation

---

## Problems Identified

### Critical Issues (Scripts Won't Run) 🔴
1. **Script 2** - References undefined functions from mrb_figure_standards.R (not sourced)
2. **Script 5** - Only conditionally sources standards file (unreliable)

### Color Inconsistencies 🟠
| Script | Current Colors | Should Be |
|--------|---------------|-----------|
| **3.abundance.R** | Hardcoded `#2b6cb0` (blue), `#e24a33` (red) | `TREATMENT_COLORS` |
| **4d.diversity.R** | `"1"="#FFD92F"` (yellow), `"3"="#8DA0CB"` (purple), `"6"="#66C2A5"` (teal) | `"1"="#E69F00"` (orange), `"3"="#56B4E9"` (sky blue), `"6"="#009E73"` (green) |
| **6.coral.R** | `"1"="#66C2A5"` (teal), `"3"="#8DA0CB"` (purple), `"6"="#FFD92F"` (yellow) - **SWAPPED!** | Correct `TREATMENT_COLORS` |
| **7.coral-caffi.R** | Same as 4d, defined **3 times** | `TREATMENT_COLORS` |

### Duplicate Function Definitions 🟡
| Script | Duplicates |
|--------|-----------|
| **4d.diversity.R** | `theme_pub()`, `save_both()`, `show_and_save()` |
| **6.coral.R** | `theme_pub()`, `save_both()`, `show_and_save()`, `save_plot()` |
| **7.coral-caffi.R** | `theme_pub()`, `save_both()`, `show_and_save()` |
| **11.nmds_permanova_cafi.R** | `save_both()`, `show_and_save()` |

### Theme Inconsistencies 🟡
- Scripts use: `theme_pub()` (local), `theme_minimal()`, `theme_bw()`, `theme_classic()`
- Should use: `theme_publication()` from mrb_figure_standards.R

### DPI Inconsistencies 🟡
- **Script 3**: Uses DPI = 300 (should be 600)
- All others: Use 600 ✅

---

## Official Standards (The Truth™)

### Color Palette (Colorblind-Friendly)
```r
TREATMENT_COLORS <- c(
  "1" = "#E69F00",   # Orange - Single coral
  "3" = "#56B4E9",   # Sky blue - Three corals
  "6" = "#009E73"    # Green - Six corals
)
```

### Theme Hierarchy
1. **theme_publication()** - Main publication theme (all scripts should use)
2. **theme_multipanel()** - For multi-panel figures
3. **theme_ordination()** - For NMDS/PCA plots
4. **theme_heatmap()** - For heatmaps

### Figure Saving
```r
save_figure(plot, filename,
           width = PUBLICATION_WIDTH_SINGLE,
           height = PUBLICATION_HEIGHT_STD,
           dpi = PUBLICATION_DPI)
```

### Scaling Functions
```r
scale_color_treatment()  # Auto-applies TREATMENT_COLORS
scale_fill_treatment()   # Auto-applies TREATMENT_COLORS
```

---

## Implementation Plan

### Phase 1: Fix Critical Issues (30 min)

#### Script 2: 2.taxonomic-coverage.R
**Action:** Add proper sourcing at top
```r
# After line 17 (after mrb_config.R source)
source("scripts/MRB/mrb_figure_standards.R")
```

#### Script 5: 5.fishes.R
**Action:** Move source to top (not inside functions)
```r
# After utils.R and mrb_config.R sources
source("scripts/MRB/mrb_figure_standards.R")
```

### Phase 2: Fix Color Inconsistencies (1 hour)

#### Script 3: 3.abundance.R
**Find & Replace:**
```r
# OLD:
col_obs = "#2b6cb0"
col_exp = "#e24a33"

# NEW:
source("scripts/MRB/mrb_figure_standards.R")
# Use scale_color_manual(values = TREATMENT_COLORS) or
# scale_color_treatment()
```

#### Script 4d: 4d.diversity.R (Line 75)
**Replace:**
```r
# OLD:
cols_trt <- c("1" = "#FFD92F", "3" = "#8DA0CB", "6" = "#66C2A5")

# NEW:
# Remove this line, use TREATMENT_COLORS from mrb_figure_standards.R
source("scripts/MRB/mrb_figure_standards.R")
```

**Then find all uses of `cols_trt` and replace with:**
```r
# OLD: scale_color_manual(values = cols_trt, ...)
# NEW: scale_color_treatment(...)

# OLD: scale_fill_manual(values = cols_trt, ...)
# NEW: scale_fill_treatment(...)
```

#### Script 6: 6.coral.R (Line 103)
**Replace:**
```r
# OLD (WRONG ORDER!):
cols_trt <- c("1" = "#66C2A5", "3" = "#8DA0CB", "6" = "#FFD92F")

# NEW:
# Remove this line, source standards file instead
source("scripts/MRB/mrb_figure_standards.R")
```

**Then replace all manual scales:**
```r
# OLD: scale_colour_manual(values = cols_trt, ...)
# NEW: scale_color_treatment(...)

# OLD: scale_fill_manual(values = cols_trt, ...)
# NEW: scale_fill_treatment(...)
```

#### Script 7: 7.coral-caffi.R (Line 631 + others)
**Replace:**
```r
# OLD (appears 3 times!):
cols_trt <- c("1" = "#FFD92F", "3" = "#8DA0CB", "6" = "#66C2A5")

# NEW:
# Remove ALL instances, source standards file once at top
source("scripts/MRB/mrb_figure_standards.R")
```

### Phase 3: Remove Duplicate Functions (1 hour)

#### Scripts: 4d, 6, 7, 11
**Remove these duplicate function definitions:**
- `theme_pub()` - Use `theme_publication()` from standards
- `save_both()` - Use `save_figure()` from standards or utils.R
- `show_and_save()` - Use `save_figure()` from standards
- `save_plot()` - Use `save_figure()` from standards

**Update function calls:**
```r
# OLD:
theme_pub()

# NEW:
theme_publication()
```

```r
# OLD:
save_both(plot, "path/file", w=8, h=6, dpi=600)
show_and_save(plot, "path/file.png", w=8, h=6)

# NEW:
save_figure(plot, "path/file",
           width = PUBLICATION_WIDTH_SINGLE,
           height = PUBLICATION_HEIGHT_STD)
# Automatically saves both .png and .pdf at 600 DPI
```

### Phase 4: Standardize Themes (30 min)

#### Update all plot theme calls:
```r
# OLD:
+ theme_minimal()
+ theme_bw()
+ theme_classic()
+ theme_pub()  # local definition

# NEW:
+ theme_publication()  # from mrb_figure_standards.R

# For specific cases:
+ theme_multipanel()   # multi-panel figures
+ theme_ordination()   # NMDS/PCA plots
+ theme_heatmap()      # heatmap figures
```

### Phase 5: Fix DPI (10 min)

#### Script 3: 3.abundance.R
**Find all ggsave() calls with dpi=300:**
```r
# OLD:
ggsave(..., dpi = 300)

# NEW:
ggsave(..., dpi = PUBLICATION_DPI)  # = 600
# Or better yet, use save_figure()
```

---

## File-by-File Checklist

### ✅ Script 1: 1.libraries.R
- No plotting code ✅

### 🔴 Script 2: 2.taxonomic-coverage.R
- [ ] Add `source("scripts/MRB/mrb_figure_standards.R")` at top
- [ ] Verify all theme_publication() calls work
- [ ] Verify all save_figure() calls work

### 🟠 Script 3: 3.abundance.R
- [ ] Add `source("scripts/MRB/mrb_figure_standards.R")` at top
- [ ] Replace hardcoded colors with TREATMENT_COLORS
- [ ] Replace theme_minimal()/theme_bw() with theme_publication()
- [ ] Change DPI from 300 to PUBLICATION_DPI (600)
- [ ] Use save_figure() instead of raw ggsave()

### 🟠 Script 4d: 4d.diversity.R
- [ ] Add `source("scripts/MRB/mrb_figure_standards.R")` at top
- [ ] Remove local `cols_trt` definition (line 75)
- [ ] Remove local `theme_pub()` definition (line 77)
- [ ] Remove local `save_both()`, `show_and_save()` functions
- [ ] Replace all `scale_*_manual(values=cols_trt)` with `scale_*_treatment()`
- [ ] Replace all `theme_pub()` with `theme_publication()`
- [ ] Update all `save_both()` calls to `save_figure()`

### 🔴 Script 5: 5.fishes.R
- [ ] Move `source("scripts/MRB/mrb_figure_standards.R")` to top (line ~17)
- [ ] Remove conditional sourcing from inside functions
- [ ] Verify all scale_*_treatment() calls work
- [ ] Verify all theme calls work

### 🟠 Script 6: 6.coral.R
- [ ] Add `source("scripts/MRB/mrb_figure_standards.R")` at top
- [ ] **Remove local `cols_trt` with WRONG colors** (line 103)
- [ ] Remove local `theme_pub()` definition (line 106)
- [ ] Remove local `save_both()`, `show_and_save()`, `save_plot()` functions
- [ ] Replace all `scale_colour_manual(values=cols_trt)` with `scale_color_treatment()`
- [ ] Replace all `scale_fill_manual(values=cols_trt)` with `scale_fill_treatment()`
- [ ] Replace all `theme_pub()` with `theme_publication()`
- [ ] Update all save functions to `save_figure()`

### 🟠 Script 7: 7.coral-caffi.R
- [ ] Add `source("scripts/MRB/mrb_figure_standards.R")` at top
- [ ] Remove ALL 3 instances of `cols_trt` definition
- [ ] Remove local `theme_pub()` definition
- [ ] Remove local `save_both()`, `show_and_save()` functions
- [ ] Replace color scales with `scale_*_treatment()`
- [ ] Replace `theme_pub()` with `theme_publication()`

### ✅ Script 8: 8.null-models.R
- Already properly sources mrb_config.R ✅
- [ ] Consider using mrb_figure_standards.R for consistency

### 🟡 Script 9: 9.coral-physio.R
- [ ] Add `source("scripts/MRB/mrb_figure_standards.R")` at top
- [ ] Replace `theme_minimal()` with `theme_publication()` or `theme_ordination()`
- [ ] Add TREATMENT_COLORS usage where applicable

### 🟡 Script 11: 11.nmds_permanova_cafi.R
- [ ] Add `source("scripts/MRB/mrb_figure_standards.R")` at top
- [ ] Remove local `save_both()`, `show_and_save()` duplicates
- [ ] Replace `theme_minimal()` with `theme_ordination()`
- [ ] Use TREATMENT_COLORS consistently

### 🟡 Script 12: 12.SLOSS.R
- [ ] Make standalone by sourcing mrb_figure_standards.R directly
- [ ] Don't rely on 4d.diversity.R for colors/themes

---

## Testing Strategy

After each script update:

1. **Syntax Check:**
   ```r
   source("scripts/MRB/[script].R", echo=TRUE)
   ```

2. **Visual Check:**
   - Verify colors match: Orange (#E69F00), Sky Blue (#56B4E9), Green (#009E73)
   - Verify all figures saved as both .png and .pdf
   - Verify DPI = 600 for all .png files

3. **Consistency Check:**
   - All plots should look similar in style
   - Fonts, margins, borders should match

---

## Expected Outcomes

### Before Standardization
- ❌ 3 different color schemes
- ❌ 5 different theme definitions
- ❌ 4 scripts with duplicate functions
- ❌ Inconsistent DPI (300 vs 600)
- ❌ 2 scripts that won't run

### After Standardization
- ✅ 1 canonical color scheme (TREATMENT_COLORS)
- ✅ Centralized theme functions (mrb_figure_standards.R)
- ✅ No duplicate function definitions
- ✅ Consistent DPI = 600
- ✅ All scripts run successfully
- ✅ Professional, publication-ready appearance

---

## Priority Order

1. **URGENT** (30 min): Fix scripts 2 & 5 (broken)
2. **HIGH** (1 hour): Fix script 6 colors (wrong mapping!)
3. **HIGH** (1 hour): Standardize scripts 3, 4d, 7 colors
4. **MEDIUM** (1 hour): Remove duplicate functions
5. **MEDIUM** (30 min): Standardize themes
6. **LOW** (10 min): Fix DPI in script 3

**Total Estimated Time: 4 hours**

---

## Next Steps

1. Review this plan
2. Start with Phase 1 (fix broken scripts)
3. Proceed systematically through phases
4. Test each script after modification
5. Commit changes with descriptive messages

---

**Note:** All changes should be backward compatible. Existing figures will look better, not different in content.
