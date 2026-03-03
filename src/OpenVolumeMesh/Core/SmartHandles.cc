#include <OpenVolumeMesh/Core/SmartHandles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.hh>

namespace OpenVolumeMesh
{

//--------------------
// Tet Smart Handles
//--------------------
SmartTetHalfFaceHandle SmartTetVertexHandle::opposite_halfface(CH ch) const {
    return kernel_.vertex_opposite_halfface(ch, *this);
}

SmartTetVertexHandle SmartTetHalfFaceHandle::opposite_vertex() const {
    return kernel_.halfface_opposite_vertex(*this);
}

//--------------------
// Hex Smart Handles
//--------------------
SmartHexHalfFaceHandle SmartHexHalfFaceHandle::opposite_halfface_in_cell() const {
    return kernel_.make_smart(kernel_.opposite_halfface_handle_in_cell(*this, cell()));
}

}
