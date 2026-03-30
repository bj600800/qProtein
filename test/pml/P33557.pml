load /Users/douzhixin/Developer/qProtein/qProtein-main/test/structure/P33557.pdb, P33557
remove solvent
color white, P33557
select cluster_0, res 96+34+69+170+172+109+79+209
show spheres, cluster_0
color wheat, cluster_0
select cluster_1, res 41+61
show spheres, cluster_1
color palegreen, cluster_1
select cluster_2, res 196+94+207
show spheres, cluster_2
color lightblue, cluster_2
select cluster_3, res 160+110
show spheres, cluster_3
color paleyellow, cluster_3
select cluster_4, res 128+137+125
show spheres, cluster_4
color lightpink, cluster_4
distance dist_salt_bridge, id 873, id 849
distance dist_salt_bridge, id 850, id 873
distance dist_salt_bridge, id 618, id 1033
distance dist_salt_bridge, id 850, id 874
distance dist_salt_bridge, id 619, id 1033
select salt_bridge, res 145+142+142+145+111+165+142+145+111+165
show sticks, salt_bridge
color atomic, salt_bridge
distance dist_disulfide_bond, id 684, id 817
select disulfide_bond, res 119+138
show sticks, disulfide_bond
color atomic, disulfide_bond
remove hydrogen
