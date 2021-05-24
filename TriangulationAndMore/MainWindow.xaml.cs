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
    ///
    

    public partial class MainWindow : Window
    {
        private PointsProcessing delaunay = new PointsProcessing();

        private GeometryModel3D SurfaceModel;
        private Model3DGroup MainModel3Dgroup = new Model3DGroup();
        // The wireframe's model.
        private GeometryModel3D VertexNormalsModel;
        private GeometryModel3D WireframeModel;

        private PerspectiveCamera TheCamera;// The change in CameraPhi when you press the up and down arrows.
        private const double CameraDPhi = 0.1;

        // The change in CameraTheta when you press the left and right arrows.
        private const double CameraDTheta = 0.05;

        // The change in CameraR when you press + or -.
        private const double CameraDR = 1;

        private double CameraPhi = Math.PI / 6.0;  // 30 degrees
        private double CameraTheta = Math.PI / 6.0;// 30 degrees
        private double CameraR = 30.0;

 
        public MainWindow()
        {
            InitializeComponent();
        }


        // Create the scene.
        // MainViewport is the Viewport3D defined
        // in the XAML code that displays everything.
        private void WindowLoaded(object sender, RoutedEventArgs e)
        {
            TheCamera = new PerspectiveCamera();
            TheCamera.FieldOfView = 60;
            MainViewport.Camera = TheCamera;
            PositionCamera();

            // Define lights.
            DefineLights();
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
            MeshGeometry3D borderMeshPlus = new MeshGeometry3D();
            MeshGeometry3D borderMeshMinus = new MeshGeometry3D();

            var maxX = 20d;
            var maxY = 40d;
            var zoom = 1d;
            var r = 7;
            var points = delaunay.GenerateGrid(maxX, maxY, zoom, r);

            

            var triangulation = delaunay.BowyerWatson(points);
            triangulation = triangulation.Where(triangle => !((triangle.Vertices[0].X == 0 && triangle.Vertices[0].Y == -(double)maxY / 4 * zoom ) || (triangle.Vertices[1].X == 0 && triangle.Vertices[1].Y == -(double)maxY / 4 * zoom ) || (triangle.Vertices[2].X == 0 && triangle.Vertices[2].Y == -(double)maxY / 4 * zoom ) ||
                                                              (triangle.Vertices[0].X == 0 && triangle.Vertices[0].Y == (double)maxY / 4 * zoom) || (triangle.Vertices[1].X == 0 && triangle.Vertices[1].Y == (double)maxY / 4 * zoom ) || (triangle.Vertices[2].X == 0 && triangle.Vertices[2].Y == (double)maxY / 4 * zoom)));


            points = points.Where(point => !((point.X == 0 && point.Y == (double)maxY / 4 * zoom) || (point.X == 0 && point.Y == -(double)maxY / 4 * zoom)));

            List<List<double>> aList = new List<List<double>>();
            List<double> bList = new List<double>();

            var plus = 5d;
            var minus = -5d;
            delaunay.GetAB(points.ToList(), out aList, out bList, plus, minus);
            var fi = delaunay.DoKachmarz(aList, bList);


            //triangulation = triangulation.Where(triangle => !((triangle.Vertices[0].X == 0 && triangle.Vertices[0].Y == 0) || (triangle.Vertices[1].X == 0 && triangle.Vertices[1].Y == 0) || (triangle.Vertices[2].X == 0 && triangle.Vertices[2].Y == 0)));

            var borderPointsMinus = points.Where(point => point.IsInnerBoundaryMinus);
            var borderPointsPlus = points.Where(point => point.IsInnerBoundaryPlus);

            //var innerPoints = points.Where(point => !point.IsBoundary);

            foreach (var point in borderPointsPlus)
            {
                Point3D p1 = new Point3D(point.X, 1, point.Y);
                int index1 = AddPoint(borderMeshPlus.Positions, p1);

                borderMeshPlus.TriangleIndices.Add(index1);
            }

            foreach (var point in borderPointsMinus)
            {
                Point3D p1 = new Point3D(point.X, 1, point.Y);
                int index1 = AddPoint(borderMeshMinus.Positions, p1);

                borderMeshMinus.TriangleIndices.Add(index1);
            }

            foreach (var triangle in triangulation)
            {
                if ((triangle.Vertices[0].X == 0 && triangle.Vertices[0].Y == -maxY / 4 * zoom) || (triangle.Vertices[1].X == 0 && triangle.Vertices[1].Y == -maxY / 4 * zoom) || (triangle.Vertices[2].X == 0 && triangle.Vertices[2].Y == -maxY / 4 * zoom) ||
                   (triangle.Vertices[0].X == 0 && triangle.Vertices[0].Y == maxY / 4 * zoom) || (triangle.Vertices[1].X == 0 && triangle.Vertices[1].Y == maxY / 4 * zoom) || (triangle.Vertices[2].X == 0 && triangle.Vertices[2].Y == maxY / 4 * zoom))
                    continue;
                var tmp = points.ToList();
                var index1 = tmp.IndexOf(tmp.FirstOrDefault(point => point.X == triangle.Vertices[0].X && point.Y == triangle.Vertices[0].Y));
                var index2 = tmp.IndexOf(tmp.FirstOrDefault(point => point.X == triangle.Vertices[1].X && point.Y == triangle.Vertices[1].Y));
                var index3 = tmp.IndexOf(tmp.FirstOrDefault(point => point.X == triangle.Vertices[2].X && point.Y == triangle.Vertices[2].Y));
                Point3D p1;
                //if (triangle.Vertices[0].IsInnerBoundaryPlus)
                //    p1 = new Point3D(triangle.Vertices[0].X, minus, triangle.Vertices[0].Y);
                //else if(triangle.Vertices[0].IsInnerBoundaryMinus)
                //     p1 = new Point3D(triangle.Vertices[0].X, plus, triangle.Vertices[0].Y);
                //else 
                //     p1 = new Point3D(triangle.Vertices[0].X, fi[index1], triangle.Vertices[0].Y);
                //Point3D p2;
                //if (triangle.Vertices[1].IsInnerBoundaryPlus)
                //    p2 = new Point3D(triangle.Vertices[1].X, minus, triangle.Vertices[1].Y);
                //else if (triangle.Vertices[1].IsInnerBoundaryMinus)
                //    p2 = new Point3D(triangle.Vertices[1].X, plus, triangle.Vertices[1].Y);
                //else
                //    p2 = new Point3D(triangle.Vertices[1].X, fi[index2], triangle.Vertices[1].Y);
                //Point3D p3;
                //if (triangle.Vertices[2].IsInnerBoundaryPlus)
                //    p3 = new Point3D(triangle.Vertices[2].X, minus, triangle.Vertices[2].Y);
                //else if (triangle.Vertices[2].IsInnerBoundaryMinus)
                //    p3 = new Point3D(triangle.Vertices[2].X, plus, triangle.Vertices[2].Y);
                //else
                //    p3 = new Point3D(triangle.Vertices[2].X, fi[index3], triangle.Vertices[2].Y);

                if (triangle.Vertices[0].IsInnerBoundaryPlus)
                    p1 = new Point3D(triangle.Vertices[0].X, plus, triangle.Vertices[0].Y);
                else if (triangle.Vertices[0].IsInnerBoundaryMinus)
                    p1 = new Point3D(triangle.Vertices[0].X, minus, triangle.Vertices[0].Y);
                else
                    p1 = new Point3D(triangle.Vertices[0].X, fi[index1], triangle.Vertices[0].Y);
                Point3D p2;
                if (triangle.Vertices[1].IsInnerBoundaryPlus)
                    p2 = new Point3D(triangle.Vertices[1].X, plus, triangle.Vertices[1].Y);
                else if (triangle.Vertices[1].IsInnerBoundaryMinus)
                    p2 = new Point3D(triangle.Vertices[1].X, minus, triangle.Vertices[1].Y);
                else
                    p2 = new Point3D(triangle.Vertices[1].X, fi[index2], triangle.Vertices[1].Y);
                Point3D p3;
                if (triangle.Vertices[2].IsInnerBoundaryPlus)
                    p3 = new Point3D(triangle.Vertices[2].X, plus, triangle.Vertices[2].Y);
                else if (triangle.Vertices[2].IsInnerBoundaryMinus)
                    p3 = new Point3D(triangle.Vertices[2].X, minus, triangle.Vertices[2].Y);
                else
                    p3 = new Point3D(triangle.Vertices[2].X, fi[index3], triangle.Vertices[2].Y);

                //Point3D p1 = new Point3D(triangle.Vertices[0].X, 0, triangle.Vertices[0].Y);
                //Point3D p2 = new Point3D(triangle.Vertices[1].X, 0, triangle.Vertices[1].Y);
                //Point3D p3 = new Point3D(triangle.Vertices[2].X, 0, triangle.Vertices[2].Y);

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
            double thickness = 0.01;
            double length = 1;

            //MeshGeometry3D vPlusNormals = borderMeshPlus.ToUpVectors(length, thickness);
            //DiffuseMaterial vnormalsPlusMaterial = new DiffuseMaterial(Brushes.Red);
            //VertexNormalsModel = new GeometryModel3D(vPlusNormals, vnormalsPlusMaterial);
            //model_group.Children.Add(VertexNormalsModel);

            //MeshGeometry3D vMinusNormals = borderMeshMinus.ToUpVectors(length, thickness);
            //DiffuseMaterial vnormalsMinusMaterial = new DiffuseMaterial(Brushes.Blue);
            //VertexNormalsModel = new GeometryModel3D(vMinusNormals, vnormalsMinusMaterial);
            //model_group.Children.Add(VertexNormalsModel);

            MeshGeometry3D wireframe = mesh.ToWireframe(thickness);
            DiffuseMaterial wireframe_material = new DiffuseMaterial(Brushes.Black);
            WireframeModel = new GeometryModel3D(wireframe, wireframe_material);
            model_group.Children.Add(WireframeModel);


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
