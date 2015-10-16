/^>/    { print }
/^[^>]/  { gsub(/[ryswmkbdhvRYSWMKBDHV]/,"n") }
/^[^>]/  { print toupper($0) }
