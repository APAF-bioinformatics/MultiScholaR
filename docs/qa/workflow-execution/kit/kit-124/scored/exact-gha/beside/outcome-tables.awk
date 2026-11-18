# markdown table: variant beside base with paired differences. awk -F'\t' -v PAIR="25 vs 7" -v SUB=full -f md.awk paired.tsv beside.tsv [paired2 beside2]
function vd(i,   s, r, v) { s = ""; for (r = 30; r <= 36; r++) { v = row[i, r]; if (v == "fail" || v == "inconclusive") s = s (s == "" ? "" : ",") (r - 29) (v == "fail" ? "F" : "i") }
    return row[i, 37] (s == "" ? "" : " (" s ")") }
function f(x, d) { if (x == "NA" || x == "") return "n/a"; return sprintf("%." d "f", x) }
FNR == 1 { next }
$6 == PAIR && $4 == SUB && NF < 30 { key = $1 "|" $2 "|" $3; P[key] = sprintf("%+.4f (%.4f)", $11, $12) " | " ($13 == "NA" ? "n/a" : sprintf("%+.4f (%.4f)", $13, $14)); next }
$6 == PAIR && $4 == SUB && NF >= 30 { key = $1 "|" $2 "|" $3; for (c = 1; c <= NF; c++) row[$7, c] = $c
    if ($7 == "base") {
        if (!hdr) { print "| Simulator, cell | Mean k (base → variant) | Type I at 0.05 | Mean FDP (SE) | Share FDP > 0.10 | Power (SE) | Verdict, base → variant (rules not passed) | Paired ΔFDP (SE) | Paired Δpower (SE) |"; print "|---|---|---|---|---|---|---|---|---|"; hdr = 1 }
        tag = $2 " " $3 ($10 == "TRUE" ? "" : ($9 == "core" ? " (void)" : " (reported)")) ($1 == "ace1f4de" ? ", kit-128 tree" : "")
        printf "| %s | %s → %s | %s → %s | %s (%s) → %s (%s) | %s → %s | %s (%s) → %s (%s) | %s → %s | %s |\n", tag,
            f(row["base", 15], 2), f(row["variant", 15], 2), f(row["base", 16], 4), f(row["variant", 16], 4),
            f(row["base", 20], 4), f(row["base", 21], 4), f(row["variant", 20], 4), f(row["variant", 21], 4),
            f(row["base", 23], 3), f(row["variant", 23], 3),
            f(row["base", 27], 3), f(row["base", 28], 4), f(row["variant", 27], 3), f(row["variant", 28], 4),
            vd("base"), vd("variant"), P[key]
    } }
