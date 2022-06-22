import graph3;
import three;
import smoothcontour3;

size(500,600,keepAspect=true);

material cubicsurf = opacity(.85) + gray;
pen bluepen = heavyblue + 1.75;

currentprojection=orthographic(34,5,-14);
//currentprojection=orthographic(-5,8,-20);
currentlight=(-1,15,22);

//has 27 real lines
real cubic(real x, real y, real z) {
  return (58313/630*x^3+13351/1260*x^2*y-2645/84*x*y^2+1997/180*y^3+1869/20*x^2*z+769/105*x*y*z-323/90*y^2*z-27/20*x*z^2-4751/90*y*z^2+12053/630*z^3+1073/120*x^2-883/120*x*y-977/252*y^2-8737/504*x*z+2171/360*y*z+541/8*z^2+1079/210*x-30977/1260*y-36767/630*z-6953/1260);
}

triple l1(real t) {return (t, - 28.4702t - 15.6456, 13.5384t + 6.99529);}
triple l2(real t) {return (t, - 10.0305t - 2.70772, - 5.43378t -.753902);}
triple l3(real t) {return (t, - 14.4001t + 4.23581, 7.09905t - 2.10375);}
triple l4(real t) {return (t, - 3.54809t + .549783, - 2.55345t+ 1.17485);}
triple l5(real t) {return (t, - 1.93403t + 1.81162, - 3.98911t + .0524936);}
triple l6(real t) {return (t, - .875448t - 1.73202, -.864684t - .266952);}
triple l7(real t) {return (t, - 4.71124t - 1.72395, 2.51475t + .716127);}
triple l8(real t) {return (t, - .992488t - .985198,-.719874t + .0735468);}
triple l9(real t) {return (t, - .751809t + 1.11422, - .929219t - 1.75255);}
triple l10(real t) {return (t, - 1.15387t -.459626, - .344849t + .534917);}
triple l11(real t) {return (t, - 1.94728t - .0431721, .808306t - .0703685);}
triple l12(real t) {return (t, - .176033t + 1.18772, -.89436t - .390832);}
triple l13(real t) {return (t, 1.67552t - .122873, - .507164t - .0414284);}
triple l14(real t) {return (t,2.25282t - .576004, - .522393t + .0984469);}
triple l15(real t) {return (t, 3.52199t - 1.60171, - .880366t + .665718);}
triple l16(real t) {return (t, - 1.781t + 4.38724, - 2.13528t + 2.08296);}
triple l17(real t) {return (t, 1.04049t + 1.98872, - .596163t + .254512);}
triple l18(real t) {return (t, 3.94658t - 1.76838, 1.54188t - .285096);}
triple l19(real t) {return (t, 3.98395t + .934384, - 1.07086t - .380084);}
triple l20(real t) {return (t, 2.72099t - .471431, 1.4078t + .529577);}
triple l21(real t) {return (t, 1.88879t + 1.51722, 2.16519t - 1.2803);}
triple l22(real t) {return (t,6.9159t + .341086, - 2.40999t - .209898);}
triple l23(real t) {return (t, 2.25034t + .12498, 1.60676t + .870069);}
triple l24(real t) {return (t,7.05555t - 4.63311, 19.2783t - 16.6281);}
triple l25(real t) {return  (t, 2.7885t + 1.91789, 1.39774t + .173715);}
triple l26(real t) {return (t,7.68822t + 5.40152, 21.0431t + 11.3636);}
triple l27(real t) {return (t, 78.7309t + 38.8121, 35.788t + 18.298);}


path3 line1 = graph(l1,-.69,-.409);
path3 line2 = graph(l2,-.668731,.128834);
path3 line3 = graph(l3,.0163759,.571928);
path3 line4 = graph(l4,-.972415,1.28232);
path3 line5 = graph(l5,-.989572,1.01589);
path3 line6 = graph(l6,-4.0347,2.59065);
path3 line7 = graph(l7,-1.21496,.483109);
path3 line8 = graph(l8,-4.02293,3.03762);
path3 line9 = graph(l9,-3.83845,2.41864);
path3 line10 = graph(l10,-3.86492,3.06826);
path3 line11 = graph(l11,-2.07632,2.03198);
path3 line12 = graph(l12,-4.00947,4.03547);
path3 line13 = graph(l13,-2.31398,2.46065);
path3 line14 = graph(l14,-1.51987,2.03124);
path3 line15 = graph(l15,-.680947,1.5905);
path3 line16 = graph(l16,.21743,2.84879);
path3 line17 = graph(l17,-4.00568,1.93301);
path3 line18 = graph(l18,-.565457,1.46161);
path3 line19 = graph(l19,-1.23857,.769492);
path3 line20 = graph(l20,-1.2968,1.64331);
path3 line21 = graph(l21,-1.25611,1.31448);
path3 line22 = graph(l22,-.627696,.529058);
path3 line23 = graph(l23,-1.83305,1.72197);
path3 line24 = graph(l24,.655045,1.07002);
path3 line25 = graph(l25,-2.12225,.746678);
path3 line26 = graph(l26,-.730103,-.349931);
path3 line27 = graph(l27,-.543777,-.442165);


//draw((0,0,0)--(1,0,0),red);
//draw((0,0,0)--(0,1,0),blue);
//draw((0,0,0)--(0,0,1),green);
draw(implicitsurface(cubic,(-4,-4,-4),(4,4,4)),surfacepen=cubicsurf);
draw(line1,bluepen);
draw(line2,bluepen);
draw(line3,bluepen);
draw(line4,bluepen);
draw(line5,bluepen);
draw(line6,bluepen);
draw(line7,bluepen);
draw(line8,bluepen);
draw(line9,bluepen);
draw(line10,bluepen);
draw(line11,bluepen);
draw(line12,bluepen);
draw(line13,bluepen);
draw(line14,bluepen);
draw(line15,bluepen);
draw(line16,bluepen);
draw(line17,bluepen);
draw(line18,bluepen);
draw(line19,bluepen);
draw(line20,bluepen);
draw(line21,bluepen);
draw(line22,bluepen);
draw(line23,bluepen);
draw(line24,bluepen);
draw(line25,bluepen);
draw(line26,bluepen);
draw(line27,bluepen);


currentprojection=
  orthographic(camera=(34.000004405770014,5.0000002824365,-14.000001205022508),
               up=(0.006688083,0.0009835416,0.01659375),
               target=(7.105427e-15,0,-5.329071e-15),
               zoom=1);
