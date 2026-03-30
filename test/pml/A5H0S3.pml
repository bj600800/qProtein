load /Users/douzhixin/Developer/qProtein/qProtein-main/test/structure/A5H0S3.pdb, A5H0S3
remove solvent
color white, A5H0S3
select cluster_0, res 56+66+44
show spheres, cluster_0
color wheat, cluster_0
select cluster_1, res 96+197+172+109+79+210+212+94
show spheres, cluster_1
color palegreen, cluster_1
select cluster_2, res 104+146
show spheres, cluster_2
color lightblue, cluster_2
select cluster_3, res 105+178+188
show spheres, cluster_3
color paleyellow, cluster_3
select cluster_4, res 110+159
show spheres, cluster_4
color lightpink, cluster_4
select cluster_5, res 126+135
show spheres, cluster_5
color palecyan, cluster_5
distance dist_salt_bridge, id 1062, id 624
distance dist_salt_bridge, id 678, id 804
distance dist_salt_bridge, id 1062, id 625
distance dist_salt_bridge, id 803, id 1037
distance dist_salt_bridge, id 804, id 679
select salt_bridge, res 164+111+117+134+164+111+134+161+134+117
show sticks, salt_bridge
color atomic, salt_bridge
remove hydrogen
