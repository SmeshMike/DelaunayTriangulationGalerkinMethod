using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Media;
using Color = System.Windows.Media.Color;

namespace TriangulationAndMore
{
    public class PointsProcessing
    {
      private double MaxX { get; set; }
        private double MaxY { get; set; }
        private IEnumerable<Triangle> border;

        public IEnumerable<Point> GeneratePoints(double maxX, double maxY, int zoom)
        {
            MaxX = maxX;
            MaxY = maxY;
            
            var point0 = new Point(-MaxX * zoom/2, -MaxY * zoom/2);
            var point1 = new Point(-MaxX * zoom/2, MaxY*zoom/2);
            var point2 = new Point(MaxX * zoom/2, MaxY * zoom/2);
            var point3 = new Point(MaxX * zoom/2, -MaxY * zoom/2);
            var points = new List<Point>(); 
            var tri1 = new Triangle(point0, point1, point3);
            var tri2 = new Triangle(point1, point2, point3);
            border = new List<Triangle>() { tri1, tri2 };
            Random random = new Random();
            for (var i = -(maxX/2); i <= maxX/2; i++)
            {
                for (var j = -(maxY/2); j <= maxY/2; j++)
                {
                    var pointX = i * zoom;
                    var pointY = j * zoom;
                    points.Add(new Point(pointX, pointY));
                }
            }


            return points;
        }

        public IEnumerable<Point> GenerateGrid(double maxX, double maxY, int zoom, int r)
        {
            MaxX = maxX;
            MaxY = maxY;

            var point0 = new Point(-MaxX * zoom / 2, -MaxY * zoom / 2);
            var point1 = new Point(-MaxX * zoom / 2, MaxY * zoom / 2);
            var point2 = new Point(MaxX * zoom / 2, MaxY * zoom / 2);
            var point3 = new Point(MaxX * zoom / 2, -MaxY * zoom / 2);
            var points = new List<Point>();
            var tri1 = new Triangle(point0, point1, point3);
            var tri2 = new Triangle(point1, point2, point3);
            border = new List<Triangle>() { tri1, tri2 };
            for (var i = -(maxX / 2); i <= maxX / 2; i++)
            {
                for (var j = -(maxY / 2); j <= maxY / 2; j++)
                {
                    if (!(((i) * (i) + (j+ maxY / 4) * (j + maxY / 4) < r * r)|| ((i) * (i) + (j - maxY / 4) * (j - maxY / 4) < r * r)))
                    {
                        var pointX = i * zoom;
                        var pointY = j * zoom;
                        points.Add(new Point(pointX, pointY));
                    }
                }
            }


            return points;
        }


        public List<Triangle> triangles;
        public List<PointF> points;

        public IEnumerable<Triangle> BowyerWatson(IEnumerable<Point> points)
        {
            //var supraTriangle = GenerateSupraTriangle();
            var triangulation = new HashSet<Triangle>(border);

            foreach (var point in points)
            {
                var badTriangles = FindBadTriangles(point, triangulation);
                var polygon = FindHoleBoundaries(badTriangles);

                foreach (var triangle in badTriangles)
                {
                    foreach (var vertex in triangle.Vertices)
                    {
                        vertex.AdjacentTriangles.Remove(triangle);
                    }
                }

                triangulation.RemoveWhere(o => badTriangles.Contains(o));

                foreach (var edge in polygon.Where(possibleEdge => possibleEdge.Point1 != point && possibleEdge.Point2 != point))
                {
                    if (!((point.X == edge.Point1.X && edge.Point1.X == edge.Point2.X) || (point.Y == edge.Point1.Y && edge.Point1.Y == edge.Point2.Y)))
                    {
                        //var triangle = new Triangle(point, edge.Point1, edge.Point2);
                        var triangle = new Triangle(edge.Point1, point, edge.Point2);
                        triangulation.Add(triangle);
                    }
                }
            }

            return triangulation;
        }

        private List<Edge> FindHoleBoundaries(ISet<Triangle> badTriangles)
        {
            var edges = new List<Edge>();
            foreach (var triangle in badTriangles)
            {
                edges.Add(new Edge(triangle.Vertices[0], triangle.Vertices[1]));
                edges.Add(new Edge(triangle.Vertices[1], triangle.Vertices[2]));
                edges.Add(new Edge(triangle.Vertices[2], triangle.Vertices[0]));
            }
            var grouped = edges.GroupBy(o => o);
            var boundaryEdges = edges.GroupBy(o => o).Where(o => o.Count() == 1).Select(o => o.First());
            return boundaryEdges.ToList();
        }

        private ISet<Triangle> FindBadTriangles(Point point, HashSet<Triangle> triangles)
        {
            var badTriangles = triangles.Where(o => o.IsPointInsideCircumcircle(point));
            return new HashSet<Triangle>(badTriangles);
        }

    }
}
