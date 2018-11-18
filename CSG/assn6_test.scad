$fn = 100;
difference() {
    cube([10, 5, 5], center=true);
    rotate([90,0,0]) translate([0,-5,0]) cylinder(10,2.5,2.5, center=true);
}
