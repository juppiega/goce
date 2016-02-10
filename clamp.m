function clamped = clamp(minVal, vals, maxVal)

clamped = min(maxVal, max(minVal, vals));

end