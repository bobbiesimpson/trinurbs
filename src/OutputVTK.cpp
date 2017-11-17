#include "OutputVTK.h"
#include "Forest.h"
#include "IParentSample.h"
#include "MultiscaleForest.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"
#include "vtkHexahedron.h"
#include "vtkCellArray.h"
#include "vtkDataSetMapper.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include <string>

namespace trinurbs
{
    void OutputVTK::outputForestGeometry(const Forest& f) const
    {
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1;     // number of cells in each parametric direction
        
        // create vtk data structures
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        
        uint sample_offset = 0;
        
        // loop over element
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto e = f.bezierElement(i);
            
            uint count = 0;
            
            // loop over sample points
            for(IParentSample isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                const auto samplept = isamplept.getCurrentPt();
                
                // get physical coordinate of sample point
                const Point3D phys_coord = e->eval(samplept.u, samplept.v, samplept.w);
//                std::cout << phys_coord << "\n";
                points->InsertPoint(sample_offset + count, phys_coord.data());
                ++count;
            }
            
            // now loop over cells and construct connectivity
            for(uint iw = 0; iw < ncell; ++iw)
            {
                for(uint iv = 0; iv < ncell; ++iv)
                {
                    for(uint iu = 0; iu < ncell; ++iu)
                    {
                        // create a hex
                        vtkSmartPointer<vtkCell> cell = vtkHexahedron::New();
                        
                        const uint index0 = iw * nsample * nsample + iv * nsample + iu;
                        const uint index2 = iw * nsample * nsample + (iv + 1) * nsample + iu + 1;
                        
                        const uint index4 = (iw + 1) * nsample * nsample + iv * nsample + iu;
                        const uint index6 = (iw + 1) * nsample * nsample + (iv + 1) * nsample + iu + 1;
                        
                        cell->GetPointIds()->SetId(0, sample_offset + index0);
                        cell->GetPointIds()->SetId(1, sample_offset +   index0 + 1);
                        cell->GetPointIds()->SetId(2, sample_offset + index2);
                        cell->GetPointIds()->SetId(3, sample_offset + index2 - 1);
                        cell->GetPointIds()->SetId(4, sample_offset + index4);
                        cell->GetPointIds()->SetId(5, sample_offset + index4 + 1 );
                        cell->GetPointIds()->SetId(6, sample_offset + index6);
                        cell->GetPointIds()->SetId(7, sample_offset + index6 - 1);
                        
                        grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds() );
                    }
                }
            }
            sample_offset += nsample * nsample * nsample;
        }
        grid->SetPoints(points);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_geometry.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        if(!writer->Write())
            error( "Cannot write vtk file" );
        
//        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
//        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
//        reader->SetFileName(fname.c_str());
//        reader->Update();
//        
//        vtkSmartPointer<vtkDataSetMapper> mapper =
//        vtkSmartPointer<vtkDataSetMapper>::New();
//        mapper->SetInputConnection(reader->GetOutputPort());
//        
//        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//        actor->SetMapper(mapper);
//        
//        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//        renderWindow->AddRenderer(renderer);
//        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//        renderWindowInteractor->SetRenderWindow(renderWindow);
//        
//        renderer->AddActor(actor);
//        renderer->SetBackground(.3, .6, .3); // Background color green
//        
//        renderWindow->Render();
//        renderWindowInteractor->Start();
    }
    
    void OutputVTK::outputMultiscaleForestGeometry(const MultiscaleForest& f) const
    {
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1;     // number of cells in each parametric direction
        
        // create vtk data structures
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
        
        uint sample_offset = 0;
        
        
        // loop over element
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto e = f.multiscaleBezierElement(i);
            
            uint count = 0;
            
            // loop over sample points
            for(IParentSample isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                auto samplept = isamplept.getCurrentPt();
                
                // get physical coordinate of sample point
                const Point3D phys_coord = e->eval(samplept.u, samplept.v, samplept.w);
                points->InsertNextPoint(phys_coord[0], phys_coord[1], phys_coord[2]);
                //points->InsertPoint(sample_offset + count, phys_coord.data());
                ++count;
            }
            
            // now loop over cells and construct connectivity
            for(uint iw = 0; iw < ncell; ++iw)
            {
                for(uint iv = 0; iv < ncell; ++iv)
                {
                    for(uint iu = 0; iu < ncell; ++iu)
                    {
                        // create a hex
                        vtkSmartPointer<vtkHexahedron> cell = vtkSmartPointer<vtkHexahedron>::New();
                        
                        const uint index0 = iw * nsample * nsample + iv * nsample + iu;
                        const uint index2 = iw * nsample * nsample + (iv + 1) * nsample + iu + 1;
                        
                        const uint index4 = (iw + 1) * nsample * nsample + iv * nsample + iu;
                        const uint index6 = (iw + 1) * nsample * nsample + (iv + 1) * nsample + iu + 1;
                        
                        cell->GetPointIds()->SetId(0, sample_offset + index0);
                        cell->GetPointIds()->SetId(1, sample_offset + index0 + 1);
                        cell->GetPointIds()->SetId(2, sample_offset + index2);
                        cell->GetPointIds()->SetId(3, sample_offset + index2 - 1);
                        cell->GetPointIds()->SetId(4, sample_offset + index4);
                        cell->GetPointIds()->SetId(5, sample_offset + index4 + 1 );
                        cell->GetPointIds()->SetId(6, sample_offset + index6);
                        cell->GetPointIds()->SetId(7, sample_offset + index6 - 1);
                        
                        cellarray->InsertNextCell(cell);
                    }
                }
            }
            sample_offset += nsample * nsample * nsample;
        }
        grid->SetPoints(points);
        grid->SetCells(VTK_HEXAHEDRON, cellarray);
        
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        
        const std::string fname = filename() + "_multiscale_geometry.vtu";
        
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        
        if(!writer->Write())
            error( "Cannot write vtk file" );
        
        // render to screen to check output valid
//        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
//        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
//        reader->SetFileName(fname.c_str());
//        reader->Update();
//        
//        vtkSmartPointer<vtkDataSetMapper> mapper =
//        vtkSmartPointer<vtkDataSetMapper>::New();
//        mapper->SetInputConnection(reader->GetOutputPort());
//        
//        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//        actor->SetMapper(mapper);
//        
//        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//        renderWindow->AddRenderer(renderer);
//        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//        renderWindowInteractor->SetRenderWindow(renderWindow);
//        
//        renderer->AddActor(actor);
//        renderer->SetBackground(.3, .6, .3); // Background color green
//        
//        renderWindow->Render();
//        renderWindowInteractor->Start();
    }
    
    void OutputVTK::outputNodalField(const Forest& f,
                                     const std::string& fieldname,
                                     const std::vector<double>& soln,
                                     const uint nlocaldof) const
    {
        
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1;     // number of cells in each parametric direction
        
        // vtk data structures
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        vtkSmartPointer<vtkDoubleArray> solndata = vtkDoubleArray::New();
        
        solndata->SetNumberOfComponents(nlocaldof);
        solndata->SetName(fieldname.c_str());
        
        // set comonent names
        for(uint idof = 0; idof < nlocaldof; ++idof)
        {
            const std::string s = "component " + std::to_string(idof);
            solndata->SetComponentName(idof, s.c_str());
        }
        
        uint sample_offset = 0;
        
        // loop over elements
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto e = f.bezierElement(i);
            const auto gbasisivec = e->globalBasisIVec();
            
            uint count = 0;
            
            // loop over sample points
            for(IParentSample isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                const auto samplept = isamplept.getCurrentPt();
                
                // get physical coordinate of sample point
                const Point3D phys_coord = e->eval(samplept.u, samplept.v, samplept.w);
//                std::cout << phys_coord << "\n";
                points->InsertPoint(sample_offset + count, phys_coord.data());
                
                // now interpolate soluiton
                const auto basisvec = e->basis(samplept.u,
                                               samplept.v,
                                               samplept.w);
                
                std::vector<double> val(nlocaldof, 0.0);
                
                for(uint ibasis = 0; ibasis < basisvec.size(); ++ibasis)
                    for(uint idof = 0; idof < nlocaldof; ++idof)
                        val[idof] += soln[gbasisivec[ibasis] * nlocaldof + idof] * basisvec[ibasis];
                
                for(uint idof = 0; idof < nlocaldof; ++idof)
                    solndata->InsertComponent(sample_offset + count, idof, val[idof]);
                ++count;
            }
            
            // now loop over cells and construct connectivity
            for(uint iw = 0; iw < ncell; ++iw)
            {
                for(uint iv = 0; iv < ncell; ++iv)
                {
                    for(uint iu = 0; iu < ncell; ++iu)
                    {
                        // create a hex
                        vtkSmartPointer<vtkCell> cell = vtkHexahedron::New();
                        
                        const uint index0 = iw * nsample * nsample + iv * nsample + iu;
                        const uint index2 = iw * nsample * nsample + (iv + 1) * nsample + iu + 1;
                        
                        const uint index4 = (iw + 1) * nsample * nsample + iv * nsample + iu;
                        const uint index6 = (iw + 1) * nsample * nsample + (iv + 1) * nsample + iu + 1;
                        
                        cell->GetPointIds()->SetId(0, sample_offset + index0);
                        cell->GetPointIds()->SetId(1, sample_offset +   index0 + 1);
                        cell->GetPointIds()->SetId(2, sample_offset + index2);
                        cell->GetPointIds()->SetId(3, sample_offset + index2 - 1);
                        cell->GetPointIds()->SetId(4, sample_offset + index4);
                        cell->GetPointIds()->SetId(5, sample_offset + index4 + 1 );
                        cell->GetPointIds()->SetId(6, sample_offset + index6);
                        cell->GetPointIds()->SetId(7, sample_offset + index6 - 1);
                        
                        grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds() );
                    }
                }
            }
            sample_offset += nsample * nsample * nsample;
        }

        grid->SetPoints(points);
        grid->GetPointData()->AddArray(solndata);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_soln.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        if(!writer->Write())
            error( "Cannot write vtk file" );
    }
    

    
}