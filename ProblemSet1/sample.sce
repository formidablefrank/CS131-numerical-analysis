//contour2d(1:10,1:10,rand(10,10),5,rect=[-2,-2,2,2])
// changing the format of the printing of the levels
//xset("fpf","%.2f")
//clf()
//contour2d(1:10,1:10,rand(10,10),5,rect=[-2,-2,2,2])

function z = el1(x, y)
    z = 2.2537948*x.^2 + 0.0063247*x*y+5.5221834*y.^2+ -1.2898102*x + -7.3773544*y
endfunction

// now an example with level numbers drawn in a legend
// Caution there are a number of tricks...
x = linspace(-2,2,30);
y = linspace(-2,2,30);
z = 2.2537948*x^2 + 0.0063247*x*y'+5.5221834*y'^2+ -1.2898102*x + -7.3773544*y'
//z = el1(x, x');
//z = cos(x')*cos(x);
//clf(); f=gcf();
//xset("fpf"," ")  // trick 1: this implies that the level numbers are not
                 //          drawn on the level curves
//f.color_map=jetcolormap(7);
// trick 2: to be able to put the legend on the right without
//          interfering with the level curves use rect with a xmax
//          value large enough 
contour2d(x,x,z,-.75:.25:.75,frameflag=3,rect=[-2,-2,2,2])
// trick 3: use legends (note that the more practical legend function
//          will not work as soon as one of the level is formed by 2 curves)  
//legends(string(-0.75:0.25:0.75),1:7,"lr");
xtitle("Graph of the orbit with perturbations")
