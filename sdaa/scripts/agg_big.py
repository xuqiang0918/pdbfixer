import re

rows = []
for line in open('/data/application/xuqiang/EPMolGen/sdaa/pdbfixer/logs/big_heavy_cmp.txt'):
    m = re.match(r'^(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s*$', line)
    if m and m.group(1) != 'pdbid':
        rows.append((m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4)),
                     int(m.group(5)), int(m.group(6))))

print('total cases:', len(rows))
ident = [r for r in rows if r[2] == 0 and r[3] == 0]
print('IDENTICAL(重原子+氢全部逐位同):', len(ident))
h_only = [r for r in rows if r[2] == 0 and r[3] > 0]
print('重原子全同,仅氢有差:', len(h_only))
heavy_diff = [r for r in rows if r[2] > 0]
print('重原子有差的 case:', len(heavy_diff))
for r in sorted(heavy_diff, key=lambda x: -x[2])[:25]:
    print('  %-6s n_heavy=%-5d heavy_max=%.3g  H_max=%.3g  重原子不一致数=%d'
          % (r[0], r[1], r[2], r[3], r[1] - r[4]))

hmax = [r[3] for r in rows]
print('H_max 分布: <=0.01: %d, 0.01-0.1: %d, 0.1-1: %d, >1: %d'
      % (sum(1 for v in hmax if v <= 0.01),
         sum(1 for v in hmax if 0.01 < v <= 0.1),
         sum(1 for v in hmax if 0.1 < v <= 1),
         sum(1 for v in hmax if v > 1)))
