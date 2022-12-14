using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	
	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.5f;                 // for collision, \mu_n 
	float mu_T = 0.5f;                          // for collision, \mu_n 

	Vector3 G = new Vector3(0.0f, -9.8f, 0.0f);


	// Use this for initialization
	void Start () 
	{		
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
	}
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}


		// In this function, update v and w by the impulse due to the collision with
		//a plane <P, N>
		void Collision_Impulse(Vector3 P, Vector3 N)
		{
			// convert quaternion to rotation matrix
			Matrix4x4 R = Matrix4x4.TRS(new Vector3(0, 0, 0), transform.rotation, new Vector3(1, 1, 1));
			// get position
			Vector3 x = transform.position;

			// If multiple points are inside plane, use average positions
			Vector3 total_ri = new Vector3(0, 0, 0);
			int n_collision_points = 0;

			// Collision detection
			Vector3[] vertices = GetComponent<MeshFilter>().mesh.vertices;
			for (int i = 0; i < vertices.Length; i++)
			{
				// get positions and velocities of each vertices given the current transform
				Vector3 ri = R * vertices[i];
				Vector3 xi = x + ri;
				Vector3 vi = v + Vector3.Cross(w, ri);

				// test if inside plane AND velocity inwards to the plane
				if (Vector3.Dot(xi - P, N) < 0 && Vector3.Dot(vi, N) < 0)
				{
					total_ri += ri;
					n_collision_points++;
				}
			}

			// Collision Handling
			if(n_collision_points > 0)
			{
				Vector3 ri = total_ri / n_collision_points;
				Vector3 vi = v + Vector3.Cross(w, ri);

				// prevent object from oscillating when staying at ground
				// if the vi y component is small after next update, set the restitution to zero
				if (Mathf.Abs(vi.y + 9.8f * dt) < 4.0f * dt) restitution = 0;

				// the moment of inertia is calculated in the ref frame
				// Now the object has an rotation R, so we should transform the moment of inertia
				Matrix4x4 inv_I = R * I_ref.inverse * R.transpose;

				// Tranverse and Normal components
				Vector3 v_n = N * Vector3.Dot(N, vi);
				Vector3 v_t = vi - v_n;

				// calculate new velocity after collision
				Vector3 v_n_temp = -v_n * restitution;
				float a = Mathf.Max(1 - mu_T * (1.0f + restitution) * v_n.magnitude / v_t.magnitude, 0.0f);
				Vector3 v_t_temp = a * v_t;
				Vector3 v_temp = v_n_temp + v_t_temp;

				// Get K matrix
				Matrix4x4 Rstar = Get_Cross_Matrix(ri);
				Matrix4x4 K = Rstar.transpose * inv_I * Rstar;
				K[0, 0] += 1.0f / mass;
				K[1, 1] += 1.0f / mass;
				K[2, 2] += 1.0f / mass;

				// calculate impulse
				
				Vector3 J = K.inverse * (v_temp-vi);

				// update velocity and angular velocity
				v += J / mass;
				w += (Vector3)(inv_I * Vector3.Cross(ri, J));

			}


		}


	// Update is called once per frame
	void Update () 
	{
		//Game Control
		if(Input.GetKey("r"))
		{
			transform.position = new Vector3 (0, 0.6f, 0);
			restitution = 0.5f;
			launched=false;
		}
		if(Input.GetKey("l"))
		{
			v = new Vector3 (5, 2, 0);
			launched=true;
		}

		// Part I: Update velocities

		v += G * dt;
		v *= linear_decay;
		w *= angular_decay;

		// Part II: Collision Impulse
		Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

		// Part III: Update position & orientation
		//Update linear status
		Vector3 x = transform.position;
		if (launched) x += dt * v;
		//Update angular status
		Quaternion q = transform.rotation;
		if (launched)
		{
			Quaternion wq = new Quaternion(w.x, w.y, w.z, 0);
			Quaternion temp_q = wq * q;
			q.x += 0.5f * dt * temp_q.x;
			q.y += 0.5f * dt * temp_q.y;
			q.z += 0.5f * dt * temp_q.z;
			q.w += 0.5f * dt * temp_q.w;
		}

		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}
