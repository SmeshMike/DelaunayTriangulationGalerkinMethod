using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;
using static TriangulationAndMore.PointsProcessing;

namespace TriangulationAndMore
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private GeometryModel3D SurfaceModel;
        private Model3DGroup MainModel3Dgroup = new Model3DGroup();
        // The wireframe's model.
        private GeometryModel3D VertexNormalsModel;
        private GeometryModel3D WireframeModel;
        private List<Triangle> triangle;
        private List<PointF> points;
        private int _pointsCount;

        private PerspectiveCamera TheCamera;// The change in CameraPhi when you press the up and down arrows.
        private const double CameraDPhi = 0.1;

        // The change in CameraTheta when you press the left and right arrows.
        private const double CameraDTheta = 0.1;

        // The change in CameraR when you press + or -.
        private const double CameraDR = 0.1;

        private double CameraPhi = Math.PI / 6.0;  // 30 degrees
        private double CameraTheta = Math.PI / 6.0;// 30 degrees
        private double CameraR = 190.0;

        const double xmin = 0;
        const double xmax = 50;
        const double dx = 1;
        const double zmin = 0;
        const double zmax = 50;
        const double dz = 1;

        public MainWindow()
        {
            InitializeComponent();
        }


        // Create the scene.
        // MainViewport is the Viewport3D defined
        // in the XAML code that displays everything.
        private void WindowLoaded(object sender, RoutedEventArgs e)
        {
            _pointsCount = Convert.ToInt32(pointsCount.Text);
            TheCamera = new PerspectiveCamera();
            TheCamera.FieldOfView = 60;
            MainViewport.Camera = TheCamera;
            PositionCamera();

            // Define lights.
            DefineLights();
            PointsProcessing pp = new PointsProcessing();
            triangle = pp.GenerateTriangles(xmin, xmax, zmin, zmax, _pointsCount);
            points = pp.points;
            // Create the model.
            DefineModel(MainModel3Dgroup);

            // Add the group of models to a ModelVisual3D.
            ModelVisual3D model_visual = new ModelVisual3D();
            model_visual.Content = MainModel3Dgroup;

            // Display the main visual in the viewportt.
            MainViewport.Children.Add(model_visual);
        }

        private void DefineModel(Model3DGroup model_group)
        {
            // Make a mesh to hold the surface.
            MeshGeometry3D mesh = new MeshGeometry3D();


            //for (double x = xmin; x <= xmax - dx; x += dx)
            //{
            //    for (double z = zmin; z <= zmax - dz; z += dx)
            //    {
            //        // Make points at the corners of the surface
            //        // over (x, z) - (x + dx, z + dz).
            //        Point3D p00 = new Point3D(x, 0, z);
            //        Point3D p10 = new Point3D(x + dx, 0, z);
            //        Point3D p01 = new Point3D(x, 0, z + dz);
            //        Point3D p11 = new Point3D(x + dx, 0, z + dz);

            //        // Add the triangles.
            //        AddTriangle(mesh, p00, p01, p11);
            //        AddTriangle(mesh, p00, p11, p10);
            //    }
            //}

            for (int i = 0; i < triangle.Count; i++)
            {
                Point3D p1 = new Point3D(triangle[i].point1.X, 0, triangle[i].point1.Y);
                Point3D p2 = new Point3D(triangle[i].point2.X, 0, triangle[i].point2.Y);
                Point3D p3 = new Point3D(triangle[i].point3.X, 0, triangle[i].point3.Y);


                // Add the triangles.
                AddTriangle(mesh, p1, p2, p3);
            }
            

            // Make the surface's material using a solid green brush.
                DiffuseMaterial surface_material = new DiffuseMaterial(Brushes.Orange);

            // Make the surface's model.
            SurfaceModel = new GeometryModel3D(mesh, surface_material);

            // Make the surface visible from both sides.
            SurfaceModel.BackMaterial = surface_material;

            // Add the model to the model groups.
            model_group.Children.Add(SurfaceModel);

            // Make a wireframe.
            double thickness = 0.3;
            double length = 5;

            MeshGeometry3D wireframe = mesh.ToWireframe(thickness);
            DiffuseMaterial wireframe_material = new DiffuseMaterial(Brushes.Black);
            WireframeModel = new GeometryModel3D(wireframe, wireframe_material);
            model_group.Children.Add(WireframeModel);

            MeshGeometry3D vnormals = mesh.ToVertexNormals(length, thickness);
            DiffuseMaterial vnormals_material = new DiffuseMaterial(Brushes.Green);
            VertexNormalsModel = new GeometryModel3D(vnormals, vnormals_material);
            model_group.Children.Add(VertexNormalsModel);
        }

        void Init()// Функция инцилизации компонентов
        {

        }

        private void TriangulationButtonClick(object sender, RoutedEventArgs e)
        {
            MainModel3Dgroup.Children.Clear();
            MainViewport.Children.Clear();
            _pointsCount = Convert.ToInt32(pointsCount.Text);
            MainViewport.Camera = TheCamera;
            PositionCamera();

            // Define lights.
            DefineLights();
            PointsProcessing pp = new PointsProcessing();
            triangle = pp.GenerateTriangles(xmin, xmax, zmin, zmax, _pointsCount);
            points = pp.points;
            // Create the model.
            DefineModel(MainModel3Dgroup);

            // Add the group of models to a ModelVisual3D.
            ModelVisual3D model_visual = new ModelVisual3D();
            model_visual.Content = MainModel3Dgroup;

            // Display the main visual in the viewportt.
            MainViewport.Children.Add(model_visual);
        }

        private void Window_KeyDown(object sender, KeyEventArgs e)
        {
            switch (e.Key)
            {
                case Key.W:
                    CameraPhi += CameraDPhi;
                    if (CameraPhi > Math.PI / 2.0) CameraPhi = Math.PI / 2.0;
                    break;
                case Key.S:
                    CameraPhi -= CameraDPhi;
                    if (CameraPhi < -Math.PI / 2.0) CameraPhi = -Math.PI / 2.0;
                    break;
                case Key.A:
                    CameraTheta += CameraDTheta;
                    break;
                case Key.D:
                    CameraTheta -= CameraDTheta;
                    break;
                case Key.Add:
                case Key.OemPlus:
                    CameraR -= CameraDR;
                    if (CameraR < CameraDR) CameraR = CameraDR;
                    break;
                case Key.Subtract:
                case Key.OemMinus:
                    CameraR += CameraDR;
                    break;
            }

            // Update the camera's position.
            PositionCamera();
        }

        // Position the camera.
        private void PositionCamera()
        {
            // Calculate the camera's position in Cartesian coordinates.
            double y = CameraR * Math.Sin(CameraPhi);
            double hyp = CameraR * Math.Cos(CameraPhi);
            double x = hyp * Math.Cos(CameraTheta);
            double z = hyp * Math.Sin(CameraTheta);
            TheCamera.Position = new Point3D(x, y, z);

            // Look toward the origin.
            TheCamera.LookDirection = new Vector3D(-x, -y, -z);

            // Set the Up direction.
            TheCamera.UpDirection = new Vector3D(0, 1, 0);

            // Console.WriteLine("Camera.Position: (" + x + ", " + y + ", " + z + ")");
        }


        private void AddTriangle(MeshGeometry3D mesh, Point3D point1, Point3D point2, Point3D point3)
        {

            // Get the points' indices.
            int index1 = AddPoint(mesh.Positions, point1);
            int index2 = AddPoint(mesh.Positions, point2);
            int index3 = AddPoint(mesh.Positions, point3);

            // Create the triangle.
            mesh.TriangleIndices.Add(index1);
            mesh.TriangleIndices.Add(index2);
            mesh.TriangleIndices.Add(index3);
        }

        private int AddPoint(Point3DCollection points, Point3D point)
        {
            // If the point is in the point dictionary,
            // return its saved index.
            if (PointDictionary.ContainsKey(point))
                return PointDictionary[point];

            // We didn't find the point. Create it.
            points.Add(point);
            PointDictionary.Add(point, points.Count - 1);
            return points.Count - 1;
        }

        // A dictionary to hold points for fast lookup.
        private Dictionary<Point3D, int> PointDictionary =
                new Dictionary<Point3D, int>();
        private void DefineLights()
        {
            AmbientLight ambient_light = new AmbientLight(Colors.Gray);
            DirectionalLight directional_light =
                    new DirectionalLight(Colors.Gray, new Vector3D(-1.0, -3.0, -2.0));
            MainModel3Dgroup.Children.Add(ambient_light);
            MainModel3Dgroup.Children.Add(directional_light);
        }
    }
}
