#include "merge_meshes.h"

#include <map>
#include <cinolib/octree.h>
#include <cinolib/merge_meshes_at_coincident_vertices.h>

#include <muselib/geometry/mesh.h>

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

template<class M, class V, class E, class P>
void merge_and_wrap_meshes(const MUSE::SurfaceMesh<M,V,E,P> &mesh0,
                           const MUSE::SurfaceMesh<M,V,E,P> &mesh1,
                           cinolib::Trimesh<> &mesh_merge)
{
    // 1. Prepara i bordi ordinati
    std::vector<uint> b0 = mesh0.get_ordered_boundary_vertices();
    std::vector<uint> b1 = mesh1.get_ordered_boundary_vertices();

    // 2. Inizializza la merge
    mesh_merge = mesh0;
    mesh_merge += mesh1;

    uint offset = mesh0.num_verts();

    // 3. Mappa corrispondenza tra vertici di bordo se coincidono
    const double tol = 1e-6;
    std::map<uint, uint> merge_map;  // mappa da vertice di mesh0 → vertice globale in mesh_merge

    for (uint vid1 : b1) {
        cinolib::vec3d v1 = mesh1.vert(vid1);
        for (uint vid0 : b0) {
            cinolib::vec3d v0 = mesh0.vert(vid0);
            if (v0.dist(v1) < tol) {
                uint global1 = offset + vid1;
                if (global1 >= mesh_merge.num_verts()) {
                    std::cerr << "ERROR: global1 index out of range\n";
                    continue;
                }
                mesh_merge.vert_merge(vid0, global1);
                // registra la corrispondenza
                merge_map[vid0] = vid1;
                break;
            }
        }
    }

    // 4. Costruisci le due liste ordinate di indici *globali* dei bordi “corrispondenti”
    std::vector<uint> or0, or1;
    size_t nb = b0.size();
    or0.reserve(nb);
    or1.reserve(nb);

    for (size_t i = 0; i < nb; ++i) {
        uint vid0 = b0[i];
        auto it = merge_map.find(vid0);
        if (it != merge_map.end()) {
            or0.push_back(vid0);
            uint vid1 = it->second;
            or1.push_back(offset + vid1);
        }
    }

    // 5. Cucitura laterale tra i border loops corrispondenti
    // Converti `mesh_merge` (MUSE) in Trimesh<>, o se `mesh_merge` è già Trimesh<>, usalo direttamente
    // Supponiamo che ci sia un modo per fare: Trimesh<> &tm = mesh_merge.as_trimesh();  // pseudocodice
    cinolib::Trimesh<> tm = mesh_merge;  // se compatibile o con conversione
    tm = add_lateral_polys(tm, or0, or1);

    // 6. Se vuoi, copia nuovamente in mesh_merge da tm, o imposta mesh_merge = tm
    // (dipende dalle tue strutture)
    mesh_merge = tm;
}

template<class M, class V, class E, class P>
void merge_and_wrap_meshes_old(const MUSE::SurfaceMesh<M,V,E,P> &mesh0, const MUSE::SurfaceMesh<M,V,E,P> &mesh1, MUSE::SurfaceMesh<M,V,E,P> &mesh_merge)
{
    std::cout << "=== Starting merge and wrap meshes ..." << std::endl;

    std::vector<uint> bverts0 = mesh0.get_ordered_boundary_vertices();
    std::vector<uint> bverts1 = mesh1.get_ordered_boundary_vertices();

    std::map <uint, uint> map; //mappa per le corrispondenze tra indici first&second
    for(uint i =0; i< bverts0.size(); i++)
        map.insert(std::pair<uint, uint> (bverts0.at(i), bverts1.at(i)));

    //std::sort(bverts1.begin(), bverts1.end(), std::greater<uint>());

    mesh_merge = mesh0;
    mesh_merge += mesh1;

    uint offset = mesh0.num_verts();

    size_t n0 = bverts0.size();

    for(uint id_vid0 = 0; id_vid0 < n0-1; id_vid0++)
    {
        cinolib::vec3d v0 = mesh0.vert(bverts0[id_vid0]);
        bool matched = false;

        //Tentativo di merging diretto
        for(uint vid1 : bverts1)
        {
            cinolib::vec3d v1 = mesh1.vert(vid1);

            if(v0.dist(v1) < 1e-06)
            {
                uint remid = offset + vid1;

                if(remid >= mesh_merge.num_verts())
                {
                    std::cout << FRED("ERROR: vert index not found!") << std::endl;
                    exit(1);
                }

                mesh_merge.vert_merge(bverts0[id_vid0], remid);
                matched = true;
                break;
            }
        }

        if(!matched)
        {
            std::map<uint,uint>::iterator it;

            uint curr0 = bverts0[id_vid0];
            uint prox0 = bverts0[id_vid0+1];

            it = map.find(curr0);
            uint curr1 = it->second;

            it = map.find(prox0);
            uint prox1 = it->second;

            //Aggiunta dei triangoli (con schema antiorario)
            std::vector<uint> polysA, polysB;
            polysA.push_back(curr0);
            polysA.push_back(curr1+offset);
            polysA.push_back(prox1+offset);

            polysB.push_back(curr0);
            polysB.push_back(prox1+offset);
            polysB.push_back(prox0);

            mesh_merge.poly_add(polysA);
            mesh_merge.poly_add(polysB);
        }
    }

    cinolib::vec3d v0 = mesh0.vert(bverts0[n0-1]);
    bool matched = false;

    //Tentativo di merging diretto
    for(uint vid1 : bverts1)
    {
        cinolib::vec3d v1 = mesh1.vert(vid1);

        if(v0.dist(v1) < 1e-06)
        {
            uint remid = offset + vid1;

            if(remid >= mesh_merge.num_verts())
            {
                std::cout << FRED("ERROR: vert index not found!") << std::endl;
                exit(1);
            }

            mesh_merge.vert_merge(bverts0[n0-1], remid);
            matched = true;
            break;
        }
    }

    if(!matched)
    {
        std::map<uint,uint>::iterator it;

        uint last0 = bverts0[n0-1];
        uint first0 = bverts0[0];

        it = map.find(last0);
        uint last1 = it->second;

        it = map.find(first0);
        uint first1 = it->second;

        //Aggiunta dei triangoli (con schema antiorario)
        std::vector<uint> polysA, polysB;
        polysA.push_back(last0);
        polysA.push_back(last1+offset);
        polysA.push_back(first1+offset);

        polysB.push_back(last0);
        polysB.push_back(first1+offset);
        polysB.push_back(first0);

        mesh_merge.poly_add(polysA);
        mesh_merge.poly_add(polysB);
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
