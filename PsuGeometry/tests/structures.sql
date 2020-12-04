SELECT distinct p.pdb_code FROM protein_structure_v1 p, protein_set_v1 s
WHERE p.pdb_code = s.pdb_code
AND s.set_name = '2019'
AND p.resolution >= 1.7
ORDER BY resolution ASC
