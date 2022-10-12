using UnityEngine;
using System.Collections;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.004f;
	float damping 	= 0.996f;
	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag=true;

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;

	private GameObject cube;
	private Bounds cube_bounds;
	private GameObject block;
	private Bounds block_bounds;


	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		Vector3[] X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}

		// Get object
		cube = GameObject.Find("Cube");
		cube_bounds = cube.GetComponent<MeshFilter>().mesh.bounds;
		block = GameObject.Find("Block");
		block_bounds = block.GetComponent<MeshFilter>().mesh.bounds;
	}

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}

	void Shallow_Wave(float[,] old_h, float[,] h, float[,] new_h)
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;

		//Step 1:
		//Compute new_h based on the shallow wave model.
		//Use Neumann boundary condition

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping;

				if (i != 0) new_h[i, j] += (h[i - 1, j] - h[i, j]) * rate;
				if (i != size - 1) new_h[i, j] += (h[i + 1, j] - h[i, j]) * rate;
				if (j != 0) new_h[i, j] += (h[i, j - 1] - h[i, j]) * rate;
				if (j != size - 1) new_h[i, j] += (h[i, j + 1] - h[i, j]) * rate;
			}
		}

		//Step 2: Block->Water coupling
		//for block 1, find the intersection with water and calculate low_h (the height of point where block is lower than water).
		Vector3 block_pos = block.transform.position;
		int block_li = (int)((block_pos.x + 5.0f) / 0.1) - 4;
		int block_ui = (int)((block_pos.x + 5.0f) / 0.1) + 4;
		int block_lj = (int)((block_pos.z + 5.0f) / 0.1) - 4;
		int block_uj = (int)((block_pos.z + 5.0f) / 0.1) + 4;
		for (int i = block_li; i <= block_ui; i++)
		{
			for (int j = block_lj; j <= block_uj; j++)
			{
				if (i >= 0 && j >= 0 && i < size && j < size)
				{
					Vector3 ray_pos_0 = block.transform.InverseTransformPoint(new Vector3(X[i * size + j].x, -10, X[i * size + j].z));
					Vector3 ray_pos_1 = block.transform.InverseTransformPoint(new Vector3(X[i * size + j].x, -9, X[i * size + j].z));
					Ray ray = new Ray(ray_pos_0, ray_pos_1 - ray_pos_0);
					float dist = 99999;
					block_bounds.IntersectRay(ray, out dist);
					low_h[i, j] = -10 + dist;
				}
			}
		}

		//then set up b and cg_mask for conjugate gradient.

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (low_h[i, j] > h[i, j])
				{
					b[i, j] = 0;
					cg_mask[i, j] = false;
					vh[i, j] = 0;
				}
				else
				{
					b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
					cg_mask[i, j] = true;
				}
			}
		}

		//Solve the Poisson equation to obtain vh (virtual height).
		Conjugate_Gradient(cg_mask, b, vh, block_li, block_ui, block_lj, block_uj);

		//for block 2, calculate low_h.

		Vector3 cube_pos = cube.transform.position;
		int cube_li = (int)((cube_pos.x + 5.0f) / 0.1) - 4;
		int cube_ui = (int)((cube_pos.x + 5.0f) / 0.1) + 4;
		int cube_lj = (int)((cube_pos.z + 5.0f) / 0.1) - 4;
		int cube_uj = (int)((cube_pos.z + 5.0f) / 0.1) + 4;
		for (int i = cube_li; i <= cube_ui; i++)
		{
			for (int j = cube_lj; j <= cube_uj; j++)
			{
				if (i >= 0 && j >= 0 && i < size && j < size)
				{
					Vector3 ray_pos_0 = cube.transform.InverseTransformPoint(new Vector3(X[i * size + j].x, -10, X[i * size + j].z));
					Vector3 ray_pos_1 = cube.transform.InverseTransformPoint(new Vector3(X[i * size + j].x, -9, X[i * size + j].z));
					Ray ray = new Ray(ray_pos_0, ray_pos_1 - ray_pos_0);
					float dist = 99999;
					cube_bounds.IntersectRay(ray, out dist);
					low_h[i, j] = -10 + dist;
				}
			}
		}
		//then set up b and cg_mask for conjugate gradient.
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (low_h[i, j] > h[i, j])
				{
					b[i, j] = 0;
					cg_mask[i, j] = false;
					vh[i, j] = 0;
				}
				else
				{
					b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
					cg_mask[i, j] = true;
				}
			}
		}

		//Solve the Poisson equation to obtain vh (virtual height).
		Conjugate_Gradient(cg_mask, b, vh, cube_li, cube_ui, cube_lj, cube_uj);

		//Diminish vh.
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (cg_mask[i, j])
					vh[i, j] *= gamma;
			}
		}

		//Update new_h by vh.
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (i != 0) new_h[i, j] += (vh[i - 1, j] - vh[i, j]) * rate;
				if (i != size - 1) new_h[i, j] += (vh[i + 1, j] - vh[i, j]) * rate;
				if (j != 0) new_h[i, j] += (vh[i, j - 1] - vh[i, j]) * rate;
				if (j != size - 1) new_h[i, j] += (vh[i, j + 1] - vh[i, j]) * rate;
			}
		}

		//Step 3
		//old_h <- h; h <- new_h;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				old_h[i, j] = h[i, j];
				h[i, j] = new_h[i, j];
			}
		}

		//Step 4: Water->Block coupling.
		GameObject Cube = GameObject.Find("Cube");
		Vector3 cube_p = Cube.transform.position;
		Mesh cube_mesh = Cube.GetComponent<MeshFilter>().mesh;

		int li = (int)((cube_p.x + 5.0f) * 10) - 4;
		int ui = (int)((cube_p.x + 5.0f) * 10) + 4;
		int lj = (int)((cube_p.z + 5.0f) * 10) - 4;
		int uj = (int)((cube_p.z + 5.0f) * 10) + 4;
		Bounds bounds = cube_mesh.bounds;

		float t = 0.004f;
		float mass = 10.0f;
		Vector3 force = new Vector3(0, -mass * 9.8f, 0);
		Vector3 torque = new Vector3(0, 0, 0);

		for (int i = li - 4; i <= ui + 4; i++)
		{ 
			for (int j = lj - 4; j <= uj + 4; j++)
			{
				if (i >= 0 && j >= 0 && i < size && j < size)
				{
					Vector3 p = new Vector3(i * 0.1f - size * 0.05f, -11, j * 0.1f - size * 0.05f);
					Vector3 q = new Vector3(i * 0.1f - size * 0.05f, -10, j * 0.1f - size * 0.05f);

					p = Cube.transform.InverseTransformPoint(p);
					q = Cube.transform.InverseTransformPoint(q);

					Ray ray = new Ray(p, q - p);
					float dist = 99999;
					bounds.IntersectRay(ray, out dist);

					if (vh[i, j] != 0)
					{
						Vector3 r = p + dist * (q - p) - cube_p;
						Vector3 f = new Vector3(0, vh[i, j], 0) * 4.0f;
						force += f;

						torque += Vector3.Cross(r, f);
					}
				}
			}
		}

		cube_v *= 0.99f;
		cube_w *= 0.99f;
		cube_v += force * t / mass;
		cube_p += cube_v * t;
		Cube.transform.position = cube_p;
		cube_w += torque * t / (100.0f * mass);
		Quaternion cube_q = Cube.transform.rotation;
		Quaternion wq = new Quaternion(cube_w.x, cube_w.y, cube_w.z, 0);
		Quaternion temp_q = wq * cube_q;
		cube_q.x += 0.5f * t * temp_q.x;
		cube_q.y += 0.5f * t * temp_q.y;
		cube_q.z += 0.5f * t * temp_q.z;
		cube_q.w += 0.5f * t * temp_q.w;
		Cube.transform.rotation = cube_q;
	}


	// Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

		//Load X.y into h.
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				h[i, j] = X[i * size + j].y;


		if (Input.GetKeyDown ("r")) 
		{
			// Add random water.
			if (Input.GetKeyDown("r"))
			{
				int i = (int)(Random.Range(0.0f, 1.0f) * size);
				int j = (int)(Random.Range(0.0f, 1.0f) * size);

				if (i < 1) i = 1;
				if (j < 1) j = 1;
				if (i >= size - 1) i = size - 2;
				if (j >= size - 1) j = size - 2;
				i = 40;
				j = 50;

				float v = 0.2f * Random.Range(0.5f, 1.0f) * 4;
				h[i, j] += v;
				h[i - 1, j] -= v / 4;
				h[i + 1, j] -= v / 4;
				h[i, j - 1] -= v / 4;
				h[i, j + 1] -= v / 4;
			}
		}

		for(int l=0; l<8; l++)
		{
			Shallow_Wave(old_h, h, new_h);
		}

		//Store h back into X.y and recalculate normal.

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				X[i * size + j].y = h[i, j];

		mesh.vertices = X;

		mesh.RecalculateNormals();


	}
}
