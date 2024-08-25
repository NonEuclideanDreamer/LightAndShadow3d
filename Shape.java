
//******************************************************************************
// author: Non-Euclidean Dreamer
// modeling 3d shapes with nodes and normal vector to draw them in light & shadow
//*********************************************************************************

//Location of Boundary and Normal vector

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

import javax.imageio.ImageIO;

public class Shape 
{
	public ArrayList<double[]> loc,n,
								br;//brightness
	public static int size=1;
	public static int start=0,width=2560*size, height=1440*size;//1600*2,height=width;//
	public static BufferedImage image=new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);
	public static DecimalFormat df=new DecimalFormat("0000");
	public static String name="sphere";
	public static double[][]zBuffer=new double[width][height];
	public static double counter=0;
	public static double res=1;
	public static void main(String[]args)
	{ 
		
		double r=height/2;
		Shape sphere=crateredSphere(new double[] {0,0,0},r,0);
	/*	sphere.addRing( new double[] {0,0,0}, r*3/4 , 2*r, Math.PI/2);
	 	sphere.addmoebius( new double[] {0,0,0}, 3*r, r/4,Math.PI/3, 1);
		sphere.addmoebius( new double[] {0,0,0}, 4*r, r/4,Math.PI/4, 2);
		sphere.addmoebius( new double[] {0,0,0}, 5*r, r/4,Math.PI/5, 3);
		sphere.addmoebius( new double[] {0,0,0}, 6*r, r/4,Math.PI/6, 4);
		sphere.addmoebius( new double[] {0,0,0}, 7*r, r/4,Math.PI/7, 5);
		sphere.parallelShadowLight(new double[] {1,0,0},0);
	*/	/*
		Shape sphere = new Shape(new ArrayList<double[]>(), new ArrayList<double[]>());
		sphere.addCone(new double[] {0,0,0}, new double[] {0,0.28,0.96},new double[] {1,0,0},500,400);
		sphere.parallelShadowLight(new double[] {-1,0,0},0);
		*/
	//	sphere.addRing( new double[] {0,0,0}, 400 , 490, Math.PI/4);
	//	sphere.addTorus(new double[] {0,0,0}, new double[] {0,0.28,0.96},new double[] {1,0,0},i,720-i);	

	//	sphere.addFloor(new double[] {-400,-400,-500},new double[] {Math.sqrt(0.8),0,Math.sqrt(0.2)},800,new double[] {0,1,0},1200);
	//	sphere.addmoebius( new double[] {0,0,0}, 400, 100,Math.PI/3, 3);
	
	//	double r=height/2;
	//	Shape sphere=new Shape(new ArrayList<double[]>(), new ArrayList<double[]>());//crateredSphere(new double[] {0,0,0},r,200);
	//	sphere.addTorus(new double[] {0,0,0},  new double[] {0,0.28,0.96},new double[] {1,0,0},r/16*11,r/10*3);
	//	sphere.addmoebius( new double[] {0,0,0}, 6*r, r/4,Math.PI*0, 7);
	//	sphere.parallelShadowLight(new double[] {-1,0,0},0);
			
		double beta=Math.PI/5,//axistilt
				cosbeta=Math.cos(beta),sinbeta=Math.sin(beta);
			
		for( int i=00;i<14400;i++)
		{	
			counter=i;
		/*	Shape sphere = new Shape(new ArrayList<double[]>(), new ArrayList<double[]>());
		//	sphere.addCone(new double[] {0,0,0}, new double[] {0,0.28,0.96},new double[] {1,0,0},i-720,1440-i);	
		//	sphere.addmoebius( new double[] {0,0,0}, i-720, 1700-i,Math.PI*0, 1);
			sphere.addmoebius( new double[] {0,0,0}, i-1080, 1930-i,Math.PI*0, 3);
			sphere.addmoebius( new double[] {0,0,0}, i-1440, 260,Math.PI*0, 5);
			sphere.parallelShadowLight(new double[] {1,0,0},0);
		*/
			double gamma=2*Math.PI*i/1500,
				cosgamma=Math.cos(gamma),singamma=Math.sin(gamma);
			System.out.println(i);
			double alpha=2*Math.PI/1440*i/**1300*/,alpha1=2*Math.PI*i/1400+1;//*i*i/8000/8000;
			sphere.parallelShadowLight(new double[] {cosgamma,0,singamma}, 1);

			sphere.parallelShadowLight(new double[] {Math.cos(alpha1),Math.sin(alpha1),0}, 2);
			clearBuffer();
			sphere.parallelDraw(new double[] {Math.cos(alpha)*cosbeta*800,Math.sin(alpha)*800,Math.cos(alpha)*sinbeta*800},
					new double[] {-Math.cos(alpha)*cosbeta,-Math.sin(alpha),-Math.cos(alpha)*sinbeta},
					new double[] {Math.sin(alpha)*cosbeta,-Math.cos(alpha),Math.sin(alpha)*sinbeta},1);
			
			image=new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);
			//sphere.wipe(0);
			sphere.wipe(1);
			sphere.wipe(2);
		}
	}
	
	private void addFloor(double[] lo, double[] dir1, double s, double[] dir2, double t) 
	{
		double[] normal=cross(dir1,dir2),mormal=times(normal,-1);
		for(int i=0;i<s*res;i++)
			for(int j=0;j<t*res;j++)
			{
				br.add(new double[] {0,0,0});
				loc.add(sum(lo,sum(sum(times(dir1,i/res),times(dir2,j/res)),normal)));
				n.add(normal);
				br.add(new double[] {0,0,0});
				loc.add(sum(lo,sum(sum(times(dir1,i/res),times(dir2,j/res)),mormal)));
				n.add(mormal);
			}
	}
	
	public void addCone(double[]center,double[]a,double[]b, double r, double l)
	{
		double[] normal=cross(a,b);
		double s=Math.sqrt(r*r+l*l);
		for(double i=-s;i<s;i+=1/res)
			for(double j=0;j<2*Math.PI;j+=Math.max(0.0001, s/r/Math.abs(i)/res))
		{
			br.add(new double[] {0,0,0});
			loc.add(sum(center,sum(sum(times(normal,i*l/s),times(a,i*r/s*Math.cos(j))),times(b,i*r/s*Math.sin(j)))));
			n.add(sum(sum(times(normal,-r/s),times(a,l/s*Math.cos(j))),times(b,l/s*Math.sin(j))));
			br.add(new double[] {0,0,0});
			loc.add(sum(center,sum(sum(times(normal,i*l/s),times(a,i*r/s*Math.cos(j))),times(b,i*r/s*Math.sin(j)))));
			n.add(sum(sum(times(normal,r/s),times(a,-l/s*Math.cos(j))),times(b,-l/s*Math.sin(j))));
		
		}
	}
	
	public void addTorus(double[] center,double[] a, double[] b, double R, double r)
	{
		double[] normal=cross(a,b);
		
		for(int i=0;i<2*Math.PI*r*res;i++)
		{
			double phi=i/res/r,
					sinphi=Math.sin(phi),
					cosphi=Math.cos(phi),
					rho=R+r*cosphi;
			for(int j=0;j<2*Math.PI*(rho)*res;j++)
			{
				double theta=j/rho/res,
						costheta=Math.cos(theta),
						sintheta=Math.sin(theta);
				br.add(new double[] {0,0,0});
				loc.add(sum(center,sum(sum(times(a,costheta*rho),times(b,sintheta*rho)),times(normal,r*sinphi))));
				n.add(sum(sum(times(a,costheta*cosphi),times(b,sintheta*cosphi)),times(normal,sinphi)));
			}
		}
	}

	public void addTrefoil(double[] center,double[] a, double[] b, double R, double r)
	{
		double[] normal=cross(a,b);
		
		for(int i=0;i<2*Math.PI*r*res;i++)
		{
			double phi=i/res/r,
					sinphi=Math.sin(phi),
					cosphi=Math.cos(phi),
					rho=R+r*cosphi;
			for(int j=0;j<2*Math.PI*(rho)*res;j++)
			{
				double theta=j/rho/res,
						costheta=Math.cos(theta),
						sintheta=Math.sin(theta);
				br.add(new double[] {0,0,0});
				loc.add(sum(center,sum(sum(times(a,costheta*rho),times(b,sintheta*rho)),times(normal,r*sinphi))));
				n.add(sum(sum(times(a,costheta*cosphi),times(b,sintheta*cosphi)),times(normal,sinphi)));
			}
		}
	}
	
	public Shape(ArrayList<double[]>l,ArrayList<double[]>normal)
	{
		loc=l; n=normal;
		br=new ArrayList<double[]>();
		for(double[] i: l)br.add(new double[] {0,0,0});
	}
	
	public static Shape crateredSphere(double[] center, double r,  int crater)
	{
		Random rand=new Random();
		ArrayList<double[]>l=new ArrayList<double[]>(),normal=new ArrayList<double[]>();
		int k=0;
		double[][]craters=new double[crater][4];//x,y,z,r
		for(int i=0;i<crater;i++)
		{
			double phi=rand.nextDouble()*2*Math.PI,
			psi=Math.asin(rand.nextDouble())*(rand.nextInt(2)*2-1)+Math.PI/2,
			d=r*8/9+rand.nextDouble(r/3);
			craters[i][3]=Math.abs(d-r)+rand.nextGaussian(r/10, r/20);
			craters[i][0]=d*Math.sin(phi)*Math.sin(psi);
			craters[i][1]=d*Math.cos(phi)*Math.sin(psi);
			craters[i][2]=d*Math.cos(psi);
			print(craters[i]);
		}
		for(int j=0;j<(int)(Math.PI*r*res);j++)
		for(int i=0;i<(int)(2*Math.PI*r*res*Math.sin(j/(r*res)));i++)
			
			{	
				int cr=-1;
				double ri=r,
						factor=r*res*Math.sin(j/(r*res));
				double[]v=new double[3];
				
				normal.add(new double[] {Math.sin(i/factor)*Math.sin(j/res/r), Math.cos(i/factor)*Math.sin(j/res/r), Math.cos(j/res/r)});
				
				for(int n=0;n<crater;n++)
				{
				if(norm(sum(times(normal.getLast(),-r),craters[n]))<craters[n][3])
					cr=n;
				}
				if(cr==-1)
				{
					double[]lc=new double[3];
					for(int m=0;m<3;m++)
					lc[m]=center[m]+r*normal.getLast()[m];
					
					l.add(lc);
				}
				else normal.removeLast();
				
			}
		for(k=0;k<crater;k++)
		for(int i=0;i<(int)(2*Math.PI*craters[k][3]*res);i++)
			for(int j=0;j<(int)(Math.PI*craters[k][3]*res);j++)
			{	
				int cr=k;
				double ri=craters[k][3];
				
				double[] no=(new double[] {Math.sin(i/res/ri)*Math.sin(j/res/ri), Math.cos(i/res/ri)*Math.sin(j/res/ri), Math.cos(j/res/ri)});
				double[]lc=sum(times(no,-ri),craters[k]);
				
				if(norm(lc)<r)
				{
				for(int n=0;n<crater;n++)if(n!=k)
				{
				if(norm(sum(lc,times(craters[n],-1)))<craters[n][3])
					cr=n;
				}
				
				if(cr==k)
				{
					normal.add(no);
					l.add(lc);
				}
				}
				
			}
		return new Shape(l,normal);
	}
	public static Shape sphere(double[] center, double r)
	{
		Random rand=new Random();
		ArrayList<double[]>l=new ArrayList<double[]>(),normal=new ArrayList<double[]>();
		int k=0;
		
		for(int i=0;i<(int)(2*Math.PI*r*res);i++)
			for(int j=0;j<(int)(Math.PI*r*res);j++)
			{	
				
				double[]v=new double[3];
				
				normal.get(k)[0]=Math.sin(i/res/r)*Math.sin(j/res/r);
				normal.get(k)[1]=Math.cos(i/res/r)*Math.sin(j/res/r);
				normal.get(k)[2]=Math.cos(j/res/r);
				
				
				for(int m=0;m<3;m++)
				{
					l.get(k)[m]=center[m]+r*normal.get(k)[m];
				}
				
				k++;
				
			}
		return new Shape(l,normal);
	}
	public void addRing(double[] center, double r1, double r2, double theta)
	{
		double costheta=Math.cos(theta),
				sintheta=Math.sin(theta);
		for(int i=(int) (r1*res);i<r2*res;i++)
		{
			for(double t=0;t<2*Math.PI;t+=1.0/i)
			{
				loc.add(sum(center,new double[] {i/res*costheta*Math.cos(t)+sintheta,i/res*Math.sin(t),i/res*sintheta*Math.cos(t)-costheta}));
				loc.add(sum(center,new double[] {i/res*costheta*Math.cos(t)-sintheta,i/res*Math.sin(t),i/res*sintheta*Math.cos(t)+costheta}));
				n.add(new double[] {sintheta,0,-costheta});
				n.add(new double[] {-sintheta,0,costheta});
				for(int j=0;j<2;j++)br.add(new double[] {0,0,0});
			}
			
		}
		for(int t=0;t<2*Math.PI*r1*res;t++)
			{
				loc.add(sum(center,new double[] {r1*costheta*Math.cos(t/res),r1*Math.sin(t/res),r1*sintheta*Math.cos(t/res)}));
				n.add(new double[] {-costheta,0,-sintheta});
				br.add(new double[] {0,0,0});
			}
		for(int t=0;t<2*Math.PI*r2*res;t++)
		{
			
				loc.add(sum(center,new double[] {r2*costheta*Math.cos(t/res),r2*Math.sin(t/res),r2*sintheta*Math.cos(t/res)}));	
				n.add(new double[] {costheta,0,sintheta});
				br.add(new double[] {0,0,0});
		}
	}
	public void addmoebius(double[] center, double R, double r, double theta,int halfs)
	{
		double costheta=Math.cos(theta),
				sintheta=Math.sin(theta);
		for(int i=(int) ((-r)*res);i<(r)*res;i++)
		{
			for(double t=0;t<4*Math.PI;t+=1.0/(R+r)/res)
			{
				double rho=R+(i)/res*Math.cos(t/2*halfs),
						x=rho*Math.cos(t),
						y=rho*Math.sin(t),
						z=(i/res)*Math.sin(t/2*halfs);
				loc.add(sum(center,new double[] {costheta*x+sintheta*z,y,-x*sintheta+z*costheta}));
				//loc.add(sum(center,new double[] {costheta*x+sintheta*z,y,-x*sintheta+z*costheta}));
				x=Math.cos(t)*Math.sin(t/2*halfs)*rho-i/res*halfs/2*Math.sin(t);
				y=Math.sin(t)*Math.sin(t/2*halfs)*rho+i/res*halfs/2*Math.cos(t);
				z=-Math.cos(t/2*halfs)*rho;
				double norm=Math.sqrt(x*x+y*y+z*z);
				n.add(new double[] {(costheta*x+sintheta*z)/norm,y/norm,(-sintheta*x+costheta*z)/norm});
			//	n.add(new double[] {-(costheta*x+sintheta*z),-y,-(-sintheta*x+costheta*z)});
				for(int j=0;j<1;j++)br.add(new double[] {0,0,0});
			}
			
		}
	/*	for(int t=0;t<2*Math.PI*r1*res;t++)
			{
				loc.add(sum(center,new double[] {r1*costheta*Math.cos(t/res),r1*Math.sin(t/res),r1*sintheta*Math.cos(t/res)}));
				n.add(new double[] {-costheta,0,-sintheta});
				br.add(new double[] {0,0,0});
			}
		for(int t=0;t<2*Math.PI*r2*res;t++)
		{
			
				loc.add(sum(center,new double[] {r2*costheta*Math.cos(t/res),r2*Math.sin(t/res),r2*sintheta*Math.cos(t/res)}));	
				n.add(new double[] {costheta,0,sintheta});
				br.add(new double[] {0,0,0});
		}*/
	}
	
	private static double norm(double[] v) {
		double out = 0;
		for(int i=0;i<v.length;i++)
		{
			out+=v[i]*v[i];
		}
		
		return Math.sqrt(out);
	}

	private static double[] sum(double[] v, double[] w) {
		double[] out=new double[v.length];
		for(int i=0;i<v.length;i++)
		{
			out[i]=v[i]+w[i];
		}
		return out;
	}

	private static double[] times(double[] v, double d) 
	{
		double[] out=new double[v.length];
		for(int i=0;i<v.length;i++)
		{
			out[i]=d*v[i];
		}
		return out;
	}

	//It makes sense in my head, alright?
	private static double[] mirvec(double a, double b, double phi, double psi, double ratio) {
		
		double sinpsi=Math.sin(psi),
				cospsi=Math.cos(psi),
				sina=Math.sin(a),
				cosa=Math.cos(a),
				sin2b=Math.sin(2*b),
				cos2b=Math.cos(2*b),
				cosaphi=Math.cos(a-phi),
				sinaphi=Math.sin(a-phi);

		return new double[] {sinpsi*(sina*cos2b*cosaphi+cosa*sinaphi)-sina*sin2b*cospsi,
							sinpsi*(cosa*cos2b*cosaphi-sina*sinaphi)-cosa*sin2b*cospsi,
							-sinpsi*sin2b*cosaphi-cos2b*cospsi
		};
	}

	public void parallelLight(double[]dir,int i)
	{
		for(int k=0;k<loc.size();k++)
		{
			double sc=times(dir,n.get(k));
			if(sc<0)
				br.get(k)[i]-=sc;
		}
	}
	
	public void parallelShadowLight(double[]dir,int i)
	{
		int j=0;
		while(Math.abs(dir[j])>.6)j++;
		double[]right=normComp(e(j),dir);
		System.out.println("ShadowLight");
		right=times(right,1.0/norm(right));
		double[]up=cross(right,dir);
		int factor=1,size=width*factor;
		double[][][]buffer=new double[size][size][2];
	//	int[][][]log=new int[size][size][800];
		ArrayList<ArrayList<Integer>>log=new ArrayList<ArrayList<Integer>>(size*size);
		for(int a=0;a<size*size;a++)log.add(new ArrayList<Integer>());
		//System.out.print(log.size());
	
		for(int k=0;k<loc.size();k++)
		{
			double sc=times(dir,n.get(k)),
					scn=res*Math.min(4,1/Math.abs(sc)),
					f=times(dir,loc.get(k))-2000;
			int		r=(int) (times(right,loc.get(k))*factor)+size/2,
					u=(int) (times(up,loc.get(k))*factor)+size/2;
		//	System.out.println(f+","+r+","+u);
			
				if(buffer[r][u][0]>f+scn)
				{
					if(sc<0)
						br.get(k)[i]-=sc;
					for(int m=0;m<2;m++)
					buffer[r][u][m]=f;
					
					for(int c=log.get(r*size+u).size()-1;c>-1;c--)
					{
						br.get(log.get(r*size+u).get(c))[i]=0;
						log.get(r*size+u).removeLast();
					}
					log.get(size*r+u).add(k);
				}
				else if(buffer[r][u][1]>f-scn)
				{
					if(sc<0)
					br.get(k)[i]-=sc;
					log.get(r*size+u).add(k);
					buffer[r][u][0]=Math.min(buffer[r][u][0], f);
					buffer[r][u][1]=Math.max(buffer[r][u][1], f);
				}
				
		}
	}

	private double[] e(int j) 
	{
		double[] out=new double[3];
		out[j]=1;
		return out;
	}

	private static double times(double[] v, double[] w) 
	{
		double out=0;
		for(int i=0;i<v.length;i++)
			out+=v[i]*w[i];
		return out;
	}
	public static void clearBuffer()
	{
		for(int i=0;i<width;i++)for(int j=0;j<height;j++)zBuffer[i][j]=10000;
	}
	
	public void parallelDraw(double[]pos, double[]dir, double[]right,double scale)
	{
		for(int k=0;k<loc.size();k++)
		{
			
			double[] v=subtract(loc.get(k),pos),
					norm=normComp(v,dir),
					par=paraComp(v,dir);
			double f=-times(n.get(k),dir);///Math.sqrt(times(v,v));
			if(f>0) 
			{
			double alpha=Math.acos( times(right,norm)/(Math.sqrt(times(right,right)*times(norm,norm))) ), d=Math.sqrt(times(norm,norm));
		//	print(v);print(norm);print(par);
	//		System.out.println("alpha="+alpha);
			if(times(dir,cross(norm,right))<0)alpha*=-1;
			int x=(int)(width/2+scale*Math.cos(alpha)*d),
					y=(int)(height/2+scale*Math.sin(alpha)*d);
			if(x>-1&&x<width&&y>-1&&y<height)
			{
				d=Math.sqrt(times(par,par));
				if(d<=zBuffer[x][y])
				{
					zBuffer[x][y]=d; 
					int[]c=new int[3];
					for(int i=0;i<3;i++)
					c[i]=(int)(br.get(k)[i]*255);//*f
					//System.out.println("blue="+c[2]);
					image.setRGB(x,y,new Color(c[0],c[1],c[2]).getRGB());
				}
			}
			}
		}
		File file=new File(name+df.format(counter+start)+".png");
		try {
			ImageIO.write(image, "png", file);
		}	catch (IOException e) {	System.out.println("IOException: Problems saving file "+name);	e.printStackTrace();}
	}
private static void print(double[] v) 
{
		System.out.print("{");
		for(int i=0;i<v.length;i++)System.out.print(v[i]+", ");
		System.out.println("}");
	}

private double[] cross(double[] a, double[] b) 
{
		return new double[]{a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
	}

	private double[] subtract(double[] a, double[] b) 
	{
		double[]out=a.clone();
		for(int i=0;i<a.length;i++)out[i]-=b[i];
		return out;
	}

	//*****************************************************
	// Returns the parallel component of vector
	//*****************************************************
	public static double[] paraComp(double[] vector,double[]dir)
	{
		double factor=times(vector,dir)/times(dir,dir);
		double[] parallel=dir.clone();
		for(int i=0;i<dir.length;i++)parallel[i]*=(factor);
		return parallel;
	}
	
	//*****************************************************
	// Returns the normal component of vector
	//*****************************************************
	public double[] normComp(double[] vector,double[]dir)
	{
		double[] normal=vector.clone(),para=paraComp(vector,dir);
		for(int i=0;i<normal.length;i++)
		normal[i]-=para[i];
		return normal;
	}
	
	//*******************************
	// takes away all light of comp i
	//*******************************
	public void wipe(int i)
	{
		for(int k=0;k<br.size();k++)
			br.get(k)[i]=0;
	}
}
