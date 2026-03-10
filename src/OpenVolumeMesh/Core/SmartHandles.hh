#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <optional>

namespace OpenVolumeMesh {

//----------------------------
// Forward Declarations
//----------------------------

class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
class PolyhedralMeshTopologyKernel;

template<class TK> class SVH;
template<class TK> class SEH;
template<class TK> class SHEH;
template<class TK> class SFH;
template<class TK> class SHFH;
template<class TK> class SCH;

template<class TopoKernel>
class SVH : public VH {
public:
    SVH(VH _vh, TopoKernel const* _kernel)
        : VH{_vh}
        , kernel_(_kernel)
    {}
    using EntityTag = typename VH::EntityTag;
    bool is_valid() const {return kernel_->is_valid(*this);}
    auto adjacent_vertices(int _max_laps = 1) const { return kernel_->vertex_vertices(*this, _max_laps); }
    auto outgoing_halfedges(int _max_laps = 1) const { return kernel_->outgoing_halfedges(*this, _max_laps); }
    auto incoming_halfedges(int _max_laps = 1) const { return kernel_->incoming_halfedges(*this, _max_laps); }
    auto edges(int _max_laps = 1) const { return kernel_->vertex_edges(*this, _max_laps); }
    auto halffaces(int _max_laps = 1) const { return kernel_->vertex_halffaces(*this, _max_laps); }
    auto point() const { return kernel_->vertex(*this); }
    auto faces(int _max_laps = 1) const { return kernel_->vertex_faces(*this, _max_laps); }
    auto cells(int _max_laps = 1) const { return kernel_->vertex_cells(*this, _max_laps); }
    size_t valence() const {return kernel_->valence(*this);}
    bool is_boundary() const {return kernel_->is_boundary(*this);}

protected:
    TopoKernel const* kernel_;
};

template<class TopoKernel>
class SEH : public EH {
public:
    SEH(EH _eh, TopoKernel const*_kernel)
        : EH{_eh}
        , kernel_(_kernel)
    {}
    bool is_valid() const {return kernel_->is_valid(*this);}
    auto vertices() const { return kernel_->edge_vertices(*this); }
    auto halffaces(int _max_laps = 1) const { return kernel_->edge_halffaces(*this, _max_laps); }
    auto faces(int _max_laps = 1) const { return kernel_->edge_faces(*this, _max_laps); }
    auto cells(int _max_laps = 1) const { return kernel_->edge_cells(*this, _max_laps); }
    SHEH<TopoKernel> halfedge_handle(int subidx) const {return make_smart(half(subidx), kernel_);}
    SHEH<TopoKernel> h0() const { return this->halfedge_handle(0);}
    SHEH<TopoKernel> h1() const { return this->halfedge_handle(1);}
    size_t valence() const {return kernel_->valence(*this);}
    bool is_boundary() const {return kernel_->is_boundary(*this);}
protected:
    TopoKernel const* kernel_;
};

template<class TopoKernel>
class SHEH : public HEH {
public:
    SHEH(HEH _heh, TopoKernel const* _kernel)
        : HEH{_heh}
        , kernel_(_kernel)
    {}
    bool is_valid() const {return kernel_->is_valid(*this);}
    SVH<TopoKernel> from() const { return make_smart(kernel_->from_vertex_handle(*this), kernel_); }
    SVH<TopoKernel> to() const { return make_smart(kernel_->to_vertex_handle(*this), kernel_); }
    SHEH<TopoKernel> opp() const { return make_smart(HEH::opp(), kernel_);}
    SEH<TopoKernel> edge() const { return make_smart(HEH::edge_handle(), kernel_);}
    SHEH<TopoKernel> opposite_handle() const {return this->opp();};
    SEH<TopoKernel> edge_handle() const {return this->edge();};
    auto halffaces(int _max_laps = 1) const { return kernel_->halfedge_halffaces(*this, _max_laps); }
    auto faces(int _max_laps = 1) const { return kernel_->halfedge_faces(*this, _max_laps); }
    auto cells(int _max_laps = 1) const { return kernel_->halfedge_cells(*this, _max_laps); }
    bool is_boundary() const {return kernel_->is_boundary(*this);}
protected:
    TopoKernel const* kernel_;
};

template<class TopoKernel>
class SFH : public FH {
public:
    SFH(FH _fh, TopoKernel const* _kernel)
        : FH{_fh}
        , kernel_(_kernel)
    {}
    bool is_valid() const {return kernel_->is_valid(*this);}
    auto vertices(int _max_laps = 1) const { return kernel_->face_vertices(*this, _max_laps); }
    auto halfedges(int _max_laps = 1) const { return kernel_->face_halfedges(*this, _max_laps); }
    auto edges(int _max_laps = 1) const { return kernel_->face_edges(*this, _max_laps); }
    auto cells() const { return kernel_->face_cells(*this); }
    SHFH<TopoKernel> h0() const { return make_smart(this->halfface_handle(0), kernel_);}
    SHFH<TopoKernel> h1() const { return make_smart(this->halfface_handle(1), kernel_);}
    size_t valence() const {return kernel_->valence(*this);}
    bool is_boundary() const {return kernel_->is_boundary(*this);}
protected:
    TopoKernel const* kernel_;
};

template<class TopoKernel>
class SHFH : public HFH {
public:
    SHFH(HFH _hfh, TopoKernel const* _kernel)
        : HFH{_hfh}
        , kernel_(_kernel)
    {}
    bool is_valid() const {return kernel_->is_valid(*this);}
    auto vertices(int _max_laps = 1) const { return kernel_->halfface_vertices(*this, _max_laps); }
    auto halfedges(int _max_laps = 1) const { return kernel_->halfface_halfedges(*this, _max_laps); }
    auto edges(int _max_laps = 1) const { return kernel_->halfface_edges(*this, _max_laps); }
    SHFH<TopoKernel> opp() const { return make_smart(HFH::opp(), kernel_);}
    SFH<TopoKernel> face() const { return make_smart(this->face_handle(), kernel_);}
    SCH<TopoKernel> cell() const { return make_smart(kernel_->incident_cell(*this), kernel_); }
    bool is_boundary() const {return kernel_->is_boundary(*this);}
protected:
    TopoKernel const* kernel_;
};

template<class TopoKernel>
class SCH : public CH {
public:
    SCH(CH _ch, TopoKernel const* _kernel)
        : CH{_ch}
        , kernel_(_kernel)
    {}
    bool is_valid() const {return kernel_->is_valid(*this);}
    auto vertices(int _max_laps = 1) const { return kernel_->cell_vertices(*this, _max_laps); }
    auto halfedges(int _max_laps = 1) const { return kernel_->cell_halfedges(*this, _max_laps); }
    auto edges(int _max_laps = 1) const { return kernel_->cell_edges(*this, _max_laps); }
    auto halffaces(int _max_laps = 1) const { return kernel_->cell_halffaces(*this, _max_laps); }
    auto faces(int _max_laps = 1) const { return kernel_->cell_faces(*this, _max_laps); }
    auto cells(int _max_laps = 1) const { return kernel_->cell_cells(*this, _max_laps); }
    size_t valence() const {return kernel_->valence(*this);}
    bool is_boundary() const {return kernel_->is_boundary(*this);}
    size_t n_vertices() const {return kernel_->n_vertices_in_cell(*this);}
protected:
    TopoKernel const* kernel_;
};

//--------------------------
// Tet Smart Handles
//--------------------------
class SmartTetVertexHandle;
class SmartTetHalfFaceHandle;

class SmartTetVertexHandle : public SVH<TetrahedralMeshTopologyKernel>
{
public:
    SmartTetVertexHandle(VH _vh, TetrahedralMeshTopologyKernel const*_kernel) :
        SVH<TetrahedralMeshTopologyKernel>(_vh, _kernel) {}

    SmartTetHalfFaceHandle opposite_halfface(CH ch) const;
};
class SmartTetHalfFaceHandle : public SHFH<TetrahedralMeshTopologyKernel>
{
public:
    SmartTetHalfFaceHandle(HFH _hfh, TetrahedralMeshTopologyKernel const*_kernel) :
        SHFH<TetrahedralMeshTopologyKernel>(_hfh, _kernel) {}

    SmartTetVertexHandle opposite_vertex() const;
};

//--------------------------
// Hex Smart Handles
//--------------------------
class SmartHexHalfFaceHandle : public SHFH<HexahedralMeshTopologyKernel>
{
public:
    SmartHexHalfFaceHandle(HFH _hfh, HexahedralMeshTopologyKernel const*_kernel) :
        SHFH<HexahedralMeshTopologyKernel>(_hfh, _kernel) {}

    SmartHexHalfFaceHandle opposite_halfface_in_cell() const;
};

//------------------------------
// Smart Handle Type Traits
//------------------------------

template<typename Entity, class TopoKernel, class Enable = void>
struct smart_handle {};
template<typename Entity, class TopoKernel>
using smart_handle_t = typename smart_handle<Entity,TopoKernel>::type;

template<class TopoKernel>
struct smart_handle<Entity::Vertex, TopoKernel> { using type = SVH<TopoKernel>;};
template<class TopoKernel>
struct smart_handle<Entity::Edge, TopoKernel> { using type = SEH<TopoKernel>;};
template<class TopoKernel>
struct smart_handle<Entity::HalfEdge, TopoKernel> { using type = SHEH<TopoKernel>;};
template<class TopoKernel>
struct smart_handle<Entity::Face, TopoKernel> { using type = SFH<TopoKernel>;};
template<class TopoKernel>
struct smart_handle<Entity::HalfFace, TopoKernel> { using type = SHFH<TopoKernel>;};
template<class TopoKernel>
struct smart_handle<Entity::Cell, TopoKernel> { using type = SCH<TopoKernel>;};

// How to set type for evrything deriving from Tet/Hex Kernel?

template<> struct smart_handle<Entity::Vertex, TetrahedralMeshTopologyKernel>
{using type = SmartTetVertexHandle;};
template<> struct smart_handle<Entity::HalfFace, TetrahedralMeshTopologyKernel>
{using type = SmartTetHalfFaceHandle;};

template<> struct smart_handle<Entity::HalfFace, HexahedralMeshTopologyKernel>
{using type = SmartHexHalfFaceHandle;};

template<class Handle, class TopoKernel>
auto make_smart(Handle _h, TopoKernel const* _kernel) {
    return smart_handle_t<tag_for_handle_t<Handle>,TopoKernel>{_h, _kernel};
}

template<typename MeshT, typename IterT>
class SmartIterT : public IterT {
public:
    SmartIterT(MeshT const&_mesh, IterT iter) 
        : IterT(std::forward<IterT>(iter))
        , mesh_(&_mesh)
    {}
    using iterator_category = typename IterT::iterator_category;
    using difference_type = typename IterT::difference_type;
    using value_type = decltype(
            make_smart(
                *std::declval<IterT>(),
                std::add_pointer_t<MeshT>()
                )); //smart_handle_t<typename IterT::Entity, MeshT>;
    // Proxy iterator: no stable address to expose via operator->.
    using pointer = void;
    using reference = value_type;
    using IterT::IterT;
    pointer operator->() const = delete;
#if 0
    {
        cur_smart_handle = mesh_->make_smart(IterT::operator*())
        return &cur_handle_;
    }
#endif

    value_type operator*() const {
        return make_smart(IterT::operator*(), mesh_);
    }
    SmartIterT& operator++() {
        IterT::operator++();
        return *this;
    }
    SmartIterT operator++(int) {
        SmartIterT before = *this;
        IterT::operator++(0);
        return before;
    }
    private:
    //value_type cur_smart_handle;
    MeshT const* mesh_ = nullptr;

};

template<typename MeshT, typename IterT>
class SmartSentinel {
public:
    SmartSentinel() = default;
    SmartSentinel(MeshT const& /*mesh*/, IterT end)
        : end_(std::forward<IterT>(end))
    {}

    friend bool operator==(SmartIterT<MeshT, IterT> const& it,
                           SmartSentinel const& s) {
        if (!s.end_) {
            return false;
        }
        return static_cast<IterT const&>(it) == *s.end_;
    }
    friend bool operator==(SmartSentinel const& s,
                           SmartIterT<MeshT, IterT> const& it) {
        return it == s;
    }
    friend bool operator!=(SmartIterT<MeshT, IterT> const& it,
                           SmartSentinel const& s) {
        return !(it == s);
    }
    friend bool operator!=(SmartSentinel const& s,
                           SmartIterT<MeshT, IterT> const& it) {
        return !(s == it);
    }

private:
    std::optional<IterT> end_;
};

template<typename _MeshT, typename _IterT>
class SmartRange
{
public:
    using MeshT = _MeshT;
    using IterT = _IterT;
    using SmartIter = SmartIterT<MeshT, _IterT>;
    using SmartSentinelT = SmartSentinel<MeshT, _IterT>;
#if defined(__cpp_lib_concepts)
    static_assert(std::input_or_output_iterator<SmartIter>);
    static_assert(std::sentinel_for<SmartSentinelT, SmartIter>);
#endif

    SmartRange(MeshT const &mesh, IterT begin, IterT end)
        : mesh_(&mesh)
        , begin_(mesh, begin)
        , end_(std::forward<IterT>(end))
        {}
    SmartIter begin() { return begin_; }
    SmartIter begin() const { return begin_; }
    SmartSentinelT end() { return SmartSentinelT(*mesh_, end_); }
    SmartSentinelT end() const { return SmartSentinelT(*mesh_, end_); }
private:
    MeshT const* mesh_;
    SmartIter begin_;
    IterT end_;
};

template<typename MeshT, typename IterT>
auto make_smart_range(const MeshT *mesh, IterT begin, IterT end)
-> SmartRange<const MeshT, IterT>
{
    return SmartRange<const MeshT, IterT>(*mesh, begin, end);
}



} // namespace OpenVolumeMesh
