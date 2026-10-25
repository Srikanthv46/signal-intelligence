# GPL reference sources: oracle / reading only

ExpressLRS (`elrs/`), DIY-Multiprotocol-TX-Module (`mpm/`), mLRS (`mlrs/`) and OpenLRSng (`olrs/`) are **GPL-3**. They
were fetched 2026-10-21 to research the RC link PHY (`../03_rc_links.md`).
- **Never copy them into `python/` or `c++/`:** the C++ library must stay GPL-free.
- Implement from the facts (rates, grids, sync words, CRC polynomials), and use these files only to cross-check.
- `cc2500_register_decoder.py` is the researcher's helper that turns CC2500 register values into rate / deviation /
  channel spacing.
