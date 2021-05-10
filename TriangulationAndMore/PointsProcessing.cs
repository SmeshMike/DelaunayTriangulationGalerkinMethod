using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO.Enumeration;
using System.Linq;
using System.Text;
using System.Windows.Documents;
using System.Windows.Media;
using System.Windows.Media.Media3D;
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

        public IEnumerable<Point> TestGenerateGrid(double maxX, double maxY, int zoom, int r)
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
                    if (((i) * (i) + (j) * (j) > r * r))
                    {
                        var point = new Point(i * zoom, j * zoom);

                        if (i == (maxX / 2) || i == -(maxX / 2) || j == -(maxY / 2) || j == (maxY / 2))
                            point.IsBoundary = true;
                        else if ((i) * (i) + (j) * (j) <= (r +0.2) * (r +0.2))
                            point.IsBoundary = true;
                        points.Add(point);
                    }
                }
            }

            points.Add(new Point(0, 0));


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

            points.Add(new Point(0, maxY / 4*zoom));
            points.Add(new Point(0, -maxY / 4 * zoom));

            return points;
        }

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

        public void GetB(List<Point> bounderyPoints, List<Point> innerPoints)
        {
            List<KeyValuePair<double,Point>> bkp = new List<KeyValuePair<double, Point>>();
            foreach (var bounderyPoint in bounderyPoints)
            {
                foreach (var innerPoint in innerPoints)
                {
                    var commonTriangles = GetCommonTriangles(bounderyPoint, innerPoint).ToList();
                    var borderTriangle = commonTriangles.Aggregate((i1, i2) => i1.BoundaryPointsCount > i2.BoundaryPointsCount ? i1 : i2);
                    var innerTriangle = commonTriangles.Aggregate((i1, i2) => i1.BoundaryPointsCount < i2.BoundaryPointsCount ? i1 : i2);
                    if (commonTriangles.Count > 0)
                    {
                        var lowerVertex1 = commonTriangles[0].Vertices.Where(coord => coord.X != bounderyPoint.X && coord.X != bounderyPoint.Y && coord.Y != innerPoint.X && coord.Y != innerPoint.Y).ToList().LastOrDefault();
                        var lowerVertex2 = commonTriangles[1].Vertices.Where(coord => coord.X != bounderyPoint.X && coord.X != bounderyPoint.Y && coord.Y != innerPoint.X && coord.Y != innerPoint.Y).ToList().LastOrDefault();
                        var m11 = new Point3D(lowerVertex1.X, lowerVertex1.Y, 0);
                        var m12 = new Point3D(lowerVertex2.X, lowerVertex2.Y, 0);
                        var m2 = new Point3D(innerPoint.X, innerPoint.Y, 1);
                        var m3 = new Point3D(bounderyPoint.X, bounderyPoint.Y, 1);
                        Triangle3D t1 = new Triangle3D(m11, m2, m3);
                        Triangle3D t2 = new Triangle3D(m12, m2, m3);
                        var b = (t1.A * t2.A + t1.B * t2.B) * (t1.S + t2.S);
                        bkp.Add(new KeyValuePair<double, Point>(b, innerPoint));
                    }
                }
            }
        }

        public void GetA(List<Point> innerPoints)
        {
            List<KeyValuePair<double, Point>> akp = new List<KeyValuePair<double, Point>>();
            foreach (var iPoint in innerPoints)
            {
                foreach (var jPoint in innerPoints)
                {
                    if (iPoint == jPoint)
                    {
                        double a = 0;
                        var triangles = iPoint.AdjacentTriangles.ToList();
                        foreach (var triangle in triangles)
                        {
                            var lowerVertex = triangle.Vertices.Where(coord => coord != iPoint).ToArray();
                            var m1 = new Point3D(lowerVertex[0].X, lowerVertex[0].Y, 0);
                            var m2 = new Point3D(lowerVertex[1].X, lowerVertex[1].Y, 0);
                            var m3 = new Point3D(iPoint.X, iPoint.Y, 1);
                            Triangle3D t3D = new Triangle3D(m1, m2, m3);
                            a += (t3D.A * t3D.A + t3D.B * t3D.B) * t3D.S;
                        }
                        akp.Add(new KeyValuePair<double, Point>(a, iPoint));
                    }
                    else if (GetCommonTriangles(iPoint, jPoint).ToList().Count > 0)
                    {
                        double a = 0;
                        var commonTriangles = GetCommonTriangles(iPoint, jPoint).ToList();
                        if (commonTriangles.Count > 0)
                        {
                            var lowerVertex1 = commonTriangles[0].Vertices.Where(coord => coord != iPoint && coord != jPoint).ToList().LastOrDefault();
                            var lowerVertex2 = commonTriangles[1].Vertices.Where(coord => coord != iPoint && coord != jPoint).ToList().LastOrDefault();
                            var m11 = new Point3D(lowerVertex1.X, lowerVertex1.Y, 0);
                            var m12 = new Point3D(lowerVertex2.X, lowerVertex2.Y, 0);
                            var m2 = new Point3D(iPoint.X, iPoint.Y, 1);
                            var m3 = new Point3D(jPoint.X, jPoint.Y, 1);
                            Triangle3D t1 = new Triangle3D(m11, m2, m3);
                            Triangle3D t2 = new Triangle3D(m12, m2, m3);
                            a += (t1.A * t2.A + t1.B * t2.B) * (t1.S + t2.S);

                        }
                    }
                }
            }
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

        public IEnumerable<Triangle> GetCommonTriangles(Point point1, Point point2)
        {
            return point1.AdjacentTriangles.Intersect(point2.AdjacentTriangles);
        }

    }
}
