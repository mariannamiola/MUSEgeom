#include "merge_meshes.h"

#include <cinolib/octree.h>
#include <cinolib/merge_meshes_at_coincident_vertices.h>

template<class M, class V, class E, class P>
void merge_meshes(const MUSE::SurfaceMesh<M,V,E,P> &mesh0, const MUSE::SurfaceMesh<M,V,E,P> &mesh1, MUSE::SurfaceMesh<M,V,E,P> &mesh_merge)
{
    std::vector<uint> bverts1 = mesh1.get_boundary_vertices();
    std::sort(bverts1.begin(), bverts1.end(), std::greater<uint>());

    mesh_merge = mesh0;
    mesh_merge += mesh1;

    for(uint vid1:bverts1)
    {
        cinolib::vec3d v1 = mesh1.vert(vid1);
        for(uint vid0:mesh0.get_boundary_vertices())
        {
            cinolib::vec3d v0 = mesh0.vert(vid0);
            if(v0 == v1)
            {
                uint remid=mesh0.num_verts()+vid1; //id vertex to remove in mesh_merge

                if(remid >= mesh_merge.num_verts())
                {
                    std::cout << FRED("ERROR: vert index not found!") << std::endl;
                    exit(1);
                }

                // std::cout << mesh0.num_verts() << std::endl;
                // std::cout << vid0 << ";" << std::setprecision(10) << mesh0.vert(vid0) << std::endl;

                // std::cout << mesh1.num_verts() << std::endl;
                // std::cout << vid1 << ";" << std::setprecision(10) << mesh1.vert(vid1) << std::endl;


                // std::cout << mesh_merge.num_verts() << std::endl;
                // std::cout << remid << ";" << std::setprecision(10) << mesh_merge.vert(remid) << std::endl;

                mesh_merge.vert_merge(vid0, remid);
            }
        }
    }
}


///
/// \brief force_boundary_match_with_normal_check
/// \param mesh0
/// \param mesh1
/// \param tol
///
template<class M, class V, class E, class P>
void force_boundary_match_with_normal_check(MUSE::SurfaceMesh<M,V,E,P> &mesh0, MUSE::SurfaceMesh<M,V,E,P> &mesh1, const double &tol)
{
    auto bverts0 = mesh0.get_ordered_boundary_vertices();
    auto bverts1 = mesh1.get_ordered_boundary_vertices();

    if (bverts0.size() != bverts1.size())
    {
        std::cerr << "=== ERROR: The boundary vertex counts do not match!" << std::endl;
        exit(0);
    }

    // Controlla se l'ordine dei bordi è opposto
    bool reversed = true;
    for (size_t i = 0; i < bverts0.size(); ++i)
    {
        double d_fwd = mesh0.vert(bverts0[i]).dist(mesh1.vert(bverts1[i]));
        double d_rev = mesh0.vert(bverts0[i]).dist(mesh1.vert(bverts1[bverts1.size() - 1 - i]));

        //double d_fwd = (mesh0.vert(bverts0[i]) - mesh1.vert(bverts1[i])).length();
        //double d_rev = (mesh0.vert(bverts0[i]) - mesh1.vert(bverts1[bverts1.size() - 1 - i])).length();

        if (d_fwd < tol)
        {
            reversed = false;
            break;
        }
        else if (d_rev < tol)
        {
            reversed = true;
            break;
        }
    }

    // Inverti l'ordine dei vertici se necessario
    if (reversed)
    {
        std::reverse(bverts1.begin(), bverts1.end());
        std::cout << "=== Boundary order reversed to match the main mesh." << std::endl;
    }

    // Forza i vertici del bordo di meshB a coincidere con quelli di meshA
    for (size_t i = 0; i < bverts0.size(); ++i)
    {
        mesh1.vert(bverts1[i]) = mesh0.vert(bverts0[i]);
    }

    std::cout << "=== Boundary vertices forced to match." << std::endl;

    // === CONTROLLA LE NORMALI ===
    mesh0.update_p_normals();
    mesh1.update_p_normals();

    std::vector<cinolib::vec3d> normals0 = mesh0.vector_poly_normals();
    std::vector<cinolib::vec3d> normals1 = mesh1.vector_poly_normals();

    uint pid0 = mesh0.adj_v2p(bverts0[0])[0];
    uint pid1 = mesh1.adj_v2p(bverts1[0])[0];


    cinolib::vec3d n0 = normals0.at(pid0);
    cinolib::vec3d n1 = normals1.at(pid1);

    // Confronta l’orientamento
    if (n0.dot(n1) < 0.0)
    {
        std::cout << "=== Normals are opposite. Flipping additional mesh ..." << std::endl;

        for(uint pid=0; pid<mesh1.num_polys(); pid++)
            mesh1.poly_flip_winding_order(pid);

        mesh1.update_p_normals();
    }
    else
    {
        std::cout << "=== Normals are compatible." << std::endl;
    }
}
