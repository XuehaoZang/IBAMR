// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// Tether (penalty) force functions.
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
System *x_solid_system, *u_solid_system;
void
tether_force_function(VectorValue<double>& F,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x_bndry, // x_bndry gives current   coordinates on the boundary mesh
                      const libMesh::Point& X_ref, // X_bndry(X_ref) gives reference coordinates on the boundary mesh
                      Elem* const elem,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double t,
                      void* /*ctx*/)
{
    // // tether_force_function() is called on elements of the boundary mesh.  Here
    // // we look up the element in the solid mesh that the current boundary
    // // element was extracted from.
    // const Elem* const interior_parent = elem->interior_parent();

    // Look up the velocity of the boundary mesh.
    const std::vector<double>& u_current = *var_data[0];

    // =========================================================================
    // 1. 运动学控制参数
    // =========================================================================
    const double f = 1.0;                       // 频率 (Hz)
    const double omega = 2.0 * M_PI * f;         
    
    const double phi_amp = 40.0 * M_PI / 180.0;  // Stroke 拍动角幅值 (40度)
    const double psi_0   = 90.0 * M_PI / 180.0;  // Pitch 俯仰中心角 (90度)
    const double psi_amp = 45.0 * M_PI / 180.0;  // Pitch 俯仰角幅值 (45度)
    
    const double t_ramp = 0.1;               // 软启动时间常数
    const double tau = 0.1;                      // 翅膀翻转时间比例
    const double C_val = 1.0 / (M_PI * tau);     

    // =========================================================================
    // 2. 绝对运动学公式计算
    // =========================================================================
    // 软启动
    double S = 1.0 - exp(-t / t_ramp);
    double S_dot = (1.0 / t_ramp) * exp(-t / t_ramp);

    double sin_wt = sin(omega * t);
    double cos_wt = cos(omega * t);
    
    // Stroke (绕 Z轴)
    double phi = phi_amp * S * sin_wt;
    double phi_dot = phi_amp * (S_dot * sin_wt + S * omega * cos_wt);

    // Pitch (绕 X轴)
    double tanh_C_cos = tanh(C_val * cos_wt);
    double tanh_C = tanh(C_val);
    
    // Pitch 角度 (中心在 90 度，t=0 时为 45 度)
    double psi = psi_0 - (psi_amp / tanh_C) * tanh_C_cos;
    
    // Pitch 速度 (链式法则求导，负负得正)
    double psi_dot = (psi_amp / tanh_C) * (1.0 - tanh_C_cos * tanh_C_cos) * (C_val * omega * sin_wt);

    // =========================================================================
    // 3. 3D 几何映射与速度传递
    // =========================================================================
    double initial_psi = M_PI / 4.0;
    double x_flat = X_ref(0);
    double y_flat = X_ref(1) * cos(-initial_psi) - X_ref(2) * sin(-initial_psi);
    double z_flat = X_ref(1) * sin(-initial_psi) + X_ref(2) * cos(-initial_psi);

    // 第二步：对平躺基准应用当前的绝对 Pitch (psi)
    double x1 = x_flat;
    double y1 = y_flat * cos(psi) - z_flat * sin(psi);
    double z1 = y_flat * sin(psi) + z_flat * cos(psi);
    
    double vx1 = 0.0;
    double vy1 = -z1 * psi_dot; 
    double vz1 =  y1 * psi_dot; 

    // 第三步：应用 Stroke (phi)
    double target_x = x1 * cos(phi) - y1 * sin(phi);
    double target_y = x1 * sin(phi) + y1 * cos(phi);
    double target_z = z1; 

    double target_vx = vx1 * cos(phi) - target_y * phi_dot;
    double target_vy = vy1 * cos(phi) + target_x * phi_dot;
    double target_vz = vz1;

    // =========================================================================
    // 4. 计算并施加 Penalty 力
    // =========================================================================
    F(0) = kappa_s * (target_x - x_bndry(0)) + eta_s * (target_vx - u_current[0]);
    F(1) = kappa_s * (target_y - x_bndry(1)) + eta_s * (target_vy - u_current[1]);
    F(2) = kappa_s * (target_z - x_bndry(2)) + eta_s * (target_vz - u_current[2]);
    
    return;
} // tether_force_function
} // namespace ModelData
using namespace ModelData;

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const string exodus_solid_filename = "solid_output.ex2"; // app_initializer->getExodusIIFilename();
        const string exodus_bndry_filename = "bndry_output.ex2"; // app_initializer->getExodusIIFilename();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = 0.5;
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                solid_mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(solid_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();

            if (elem_type == "TRI6") solid_mesh.all_second_order();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        // {
        //     // NOTE: number of segments along boundary is 4*2^r.
        //     const double num_circum_segments = 2.0 * M_PI * R / ds;
        //     const int r = log2(0.25 * num_circum_segments);
        //     MeshTools::Generation::build_sphere(solid_mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        // }

        // // Ensure nodes on the surface are on the analytic boundary.
        // MeshBase::element_iterator el_end = solid_mesh.elements_end();
        // for (MeshBase::element_iterator el = solid_mesh.elements_begin(); el != el_end; ++el)
        // {
        //     Elem* const elem = *el;
        //     for (unsigned int side = 0; side < elem->n_sides(); ++side)
        //     {
        //         const bool at_mesh_bdry = !elem->neighbor_ptr(side);
        //         if (!at_mesh_bdry) continue;
        //         for (unsigned int k = 0; k < elem->n_nodes(); ++k)
        //         {
        //             if (!elem->is_node_on_side(k, side)) continue;
        //             Node& n = elem->node_ref(k);
        //             n = R * n.unit();
        //         }
        //     }
        // }
        // solid_mesh.prepare_for_use();

        // MeshTools::Generation::build_cube(solid_mesh, 16, 8, 2, 
        //                                   0.0, 1.0,      // X: Spanwise (展向)
        //                                   -0.25, 0.25,     // Y: Chordwise (弦向)
        //                                   -0.0625, 0.0625,   // Z: Thickness (厚度)
        //                                   HEX8);

        // 按照真实果蝇机翼比例 (250:75:7) 生成代理模型
        // 展长 1.0, 弦长 0.3, 厚度 0.028
        // 固体网格间距 ds 严格控制在 ~0.03 左右，完美匹配 0.0625 的流体网格
        MeshTools::Generation::build_cube(solid_mesh, 
                                          32, 10, 1,             // 切分数量：X(32段), Y(10段), Z(1段-上下两层皮)
                                          0.0, 1.0,              // X: Spanwise (展向)
                                          -0.15, 0.15,           // Y: Chordwise (弦向)
                                          -0.014, 0.014,         // Z: Thickness (厚度)
                                          HEX8);
        // solid_mesh.read("wingSolid.msh");

        // Pre-pitch 45 deg
        const double initial_psi = M_PI / 4.0;
        for (MeshBase::node_iterator it = solid_mesh.nodes_begin(); it != solid_mesh.nodes_end(); ++it)
        {
            Node* n = *it;
            double y_orig = (*n)(1);
            double z_orig = (*n)(2);
            // 绕 X 轴旋转矩阵
            (*n)(1) = y_orig * cos(initial_psi) - z_orig * sin(initial_psi);
            (*n)(2) = y_orig * sin(initial_psi) + z_orig * cos(initial_psi);
        }

        // 准备使用网格
        solid_mesh.prepare_for_use();

        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &solid_mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1, SystemData(ib_method_ops->getVelocitySystemName(), vars));
        IBFEMethod::LagBodyForceFcnData body_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagBodyForceFunction(body_fcn_data);
        EquationSystems* bndry_equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        // Setup solid systems.
        std::unique_ptr<EquationSystems> solid_equation_systems = std::make_unique<EquationSystems>(solid_mesh);
        x_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("position");
        u_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("velocity");
        Order order = FIRST;
        FEFamily family = LAGRANGE;
        for (int d = 0; d < NDIM; ++d)
        {
            x_solid_system->add_variable("X_" + std::to_string(d), order, family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            u_solid_system->add_variable("U_" + std::to_string(d), order, family);
        }
        solid_equation_systems->init();

        // Set up the position vector.
        {
            MeshBase& mesh = solid_equation_systems->get_mesh();
            System& X_system = solid_equation_systems->get_system("position");
            const unsigned int X_sys_num = X_system.number();
            NumericVector<double>& X_coords = *X_system.solution;
            for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
            {
                Node* n = *it;
                if (n->n_vars(X_sys_num))
                {
                    TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
                    const libMesh::Point& X = *n;
                    libMesh::Point x = X;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int dof_index = n->dof_number(X_sys_num, d, 0);
                        X_coords.set(dof_index, x(d));
                    }
                }
            }
            X_coords.close();
            X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
            copy_and_synch(*X_system.solution, *X_system.current_local_solution);
        }

        x_solid_system->assemble_before_solve = false;
        x_solid_system->assemble();

        u_solid_system->assemble_before_solve = false;
        u_solid_system->assemble();

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_solid_io =
            uses_exodus ? std::make_unique<ExodusII_IO>(solid_mesh) : nullptr;
        std::unique_ptr<ExodusII_IO> exodus_bndry_io =
            uses_exodus ? std::make_unique<ExodusII_IO>(boundary_mesh) : nullptr;

        // Check to see if this is a restarted run to append current exodus files
        if (uses_exodus)
        {
            const bool from_restart = RestartManager::getManager()->isFromRestart();
            exodus_solid_io->append(from_restart);
            exodus_bndry_io->append(from_restart);
        }

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_solid_io->write_timestep(
                    exodus_solid_filename, *solid_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_bndry_io->write_timestep(
                    exodus_bndry_filename, *bndry_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // ==============================================================
        // 创建用于记录受力的 CSV 文件 (包含 Manual 和 System 对比)
        // ==============================================================
        std::ofstream force_file;
        if (SAMRAI_MPI::getRank() == 0) {
            force_file.open("aero_forces_compare.csv");
            // 写入表头: 手动积分的力 (man) 和 系统提取的力 (sys)
            force_file << "Time,Fx_man,Fy_man,Fz_man,Fx_sys,Fy_sys,Fz_sys\n"; 
        }
        
        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            // Setup the position and velocity vector. --> update VisIt viz
            {
                MeshBase& mesh = solid_equation_systems->get_mesh();
                System& X_system = solid_equation_systems->get_system("position");
                const unsigned int X_sys_num = X_system.number();
                NumericVector<double>& X_coords = *X_system.solution;

                // System& U_system = solid_equation_systems->get_system("velocity");
                // const unsigned int U_sys_num = U_system.number();
                // NumericVector<double>& U_coords = *U_system.solution;
                
                const double f = 1.0;
                const double omega = 2.0 * M_PI * f;
                const double phi_amp = 40.0 * M_PI / 180.0;
                const double psi_0   = 90.0 * M_PI / 180.0;
                const double psi_amp = 45.0 * M_PI / 180.0;
                const double t_ramp = 0.1;
                const double tau = 0.1;
                const double C_val = 1.0 / (M_PI * tau);

                double t = loop_time;
                double S = 1.0 - exp(-t / t_ramp);
                double phi = phi_amp * S * sin(omega * t);
                double tanh_C_cos = tanh(C_val * cos(omega * t));
                double psi = psi_0 - (psi_amp / tanh(C_val)) * tanh_C_cos;

                double initial_psi = M_PI / 4.0; // 解旋角

                for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
                {
                    Node* n = *it;
                    if (n->n_vars(X_sys_num))
                    {
                        const libMesh::Point& X_ref = *n; 
                        // 解旋
                        double x_flat = X_ref(0);
                        double y_flat = X_ref(1) * cos(-initial_psi) - X_ref(2) * sin(-initial_psi);
                        double z_flat = X_ref(1) * sin(-initial_psi) + X_ref(2) * cos(-initial_psi);

                        // 绝对 Pitch
                        double x1 = x_flat;
                        double y1 = y_flat * cos(psi) - z_flat * sin(psi);
                        double z1 = y_flat * sin(psi) + z_flat * cos(psi);

                        // 绝对 Stroke
                        double target_x = x1 * cos(phi) - y1 * sin(phi);
                        double target_y = x1 * sin(phi) + y1 * cos(phi);
                        double target_z = z1;

                        const int dof_x = n->dof_number(X_sys_num, 0, 0);
                        const int dof_y = n->dof_number(X_sys_num, 1, 0);
                        const int dof_z = n->dof_number(X_sys_num, 2, 0);

                        X_coords.set(dof_x, target_x);
                        X_coords.set(dof_y, target_y);
                        X_coords.set(dof_z, target_z);
                    }
                }

                X_coords.close();
                X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
                copy_and_synch(X_coords, *X_system.current_local_solution);
                // U_coords.close();
                // U_system.get_dof_map().enforce_constraints_exactly(U_system, &U_coords);
                // copy_and_synch(U_coords, *U_system.current_local_solution);
            }

            // ==============================================================
            // 双通道计算/提取机翼受到的真实空气动力 (升力与阻力)
            // ==============================================================
            {
                MeshBase& mesh = solid_equation_systems->get_mesh();
                
                // 1. 获取物理位移系统
                System& coord_system = bndry_equation_systems->get_system("IB coordinates system");
                NumericVector<double>& actual_coords = *coord_system.current_local_solution;
                const unsigned int coord_sys_num = coord_system.number();

                // 2. 获取物理受力系统 (这是你刚才用探针查出来的宝藏！)
                System& force_system = bndry_equation_systems->get_system("IB force system");
                NumericVector<double>& actual_forces = *force_system.current_local_solution;
                const unsigned int force_sys_num = force_system.number();

                // 运动学参数 (保持与 tether 里面一致)
                const double f = 1.0;
                const double omega = 2.0 * M_PI * f;
                const double phi_amp = 40.0 * M_PI / 180.0;
                const double psi_0   = 90.0 * M_PI / 180.0;
                const double psi_amp = 45.0 * M_PI / 180.0;
                const double t_ramp = 0.1;
                const double tau = 0.1;
                const double C_val = 1.0 / (M_PI * tau);
                const double initial_psi = M_PI / 4.0; 
                
                const double KAPPA_S = 1.0e6; // 仅用于 Manual 计算核对

                double t = loop_time;
                double S = 1.0 - exp(-t / t_ramp);
                double phi = phi_amp * S * sin(omega * t);
                double tanh_C_cos = tanh(C_val * cos(omega * t));
                double psi = psi_0 - (psi_amp / tanh(C_val)) * tanh_C_cos;

                // 两个通道的受力累加器
                double local_Fx_man = 0.0, local_Fy_man = 0.0, local_Fz_man = 0.0;
                double local_Fx_sys = 0.0, local_Fy_sys = 0.0, local_Fz_sys = 0.0;

                for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
                {
                    Node* n = *it;
                    // 确保节点同时在这两个系统里有注册
                    if (n->n_vars(coord_sys_num) && n->n_vars(force_sys_num))
                    {
                        // ----- 通道 A: 手动算位移差 (Manual) -----
                        const libMesh::Point& X_ref = *n; 
                        double x_flat = X_ref(0);
                        double y_flat = X_ref(1) * cos(-initial_psi) - X_ref(2) * sin(-initial_psi);
                        double z_flat = X_ref(1) * sin(-initial_psi) + X_ref(2) * cos(-initial_psi);

                        double x1 = x_flat;
                        double y1 = y_flat * cos(psi) - z_flat * sin(psi);
                        double z1 = y_flat * sin(psi) + z_flat * cos(psi);
                        
                        double target_x = x1 * cos(phi) - y1 * sin(phi);
                        double target_y = x1 * sin(phi) + y1 * cos(phi);
                        double target_z = z1;

                        const int dof_cx = n->dof_number(coord_sys_num, 0, 0);
                        const int dof_cy = n->dof_number(coord_sys_num, 1, 0);
                        const int dof_cz = n->dof_number(coord_sys_num, 2, 0);
                        double actual_x = actual_coords(dof_cx);
                        double actual_y = actual_coords(dof_cy);
                        double actual_z = actual_coords(dof_cz);

                        local_Fx_man += -KAPPA_S * (target_x - actual_x);
                        local_Fy_man += -KAPPA_S * (target_y - actual_y);
                        local_Fz_man += -KAPPA_S * (target_z - actual_z);

                        // ----- 通道 B: 直接读取底层体力密度 (System) -----
                        const int dof_fx = n->dof_number(force_sys_num, 0, 0);
                        const int dof_fy = n->dof_number(force_sys_num, 1, 0);
                        const int dof_fz = n->dof_number(force_sys_num, 2, 0);
                        
                        // IBAMR 内部存储的是“固体施加给流体的力”。
                        // 根据牛顿第三定律，机翼受到的气动力是它的相反数，所以加负号。
                        local_Fx_sys += -actual_forces(dof_fx);
                        local_Fy_sys += -actual_forces(dof_fy);
                        local_Fz_sys += -actual_forces(dof_fz);
                    }
                }

                // 乘以节点控制体积 (黎曼和求积分)
                // Span(1.0) * Chord(0.3) * Thickness(0.028)
                double V_total = 1.0 * 0.3 * 0.028; 
                double total_nodes = mesh.n_nodes(); 
                double dV = V_total / total_nodes;

                local_Fx_man *= dV; local_Fy_man *= dV; local_Fz_man *= dV;
                local_Fx_sys *= dV; local_Fy_sys *= dV; local_Fz_sys *= dV;

                // MPI 并行规约合并 (6 个变量一起合并)
                double global_F[6] = {local_Fx_man, local_Fy_man, local_Fz_man, 
                                      local_Fx_sys, local_Fy_sys, local_Fz_sys};
                SAMRAI_MPI::sumReduction(global_F, 6);

                // 写入 CSV
                if (SAMRAI_MPI::getRank() == 0) {
                    force_file << loop_time << "," 
                               << global_F[0] << "," << global_F[1] << "," << global_F[2] << ","
                               << global_F[3] << "," << global_F[4] << "," << global_F[5] << "\n";
                    force_file.flush(); 
                }
            }

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_solid_io->write_timestep(exodus_solid_filename,
                                                    *solid_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                    exodus_bndry_io->write_timestep(exodus_bndry_filename,
                                                    *bndry_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main
