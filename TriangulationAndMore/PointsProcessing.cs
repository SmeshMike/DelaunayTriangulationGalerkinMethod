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

        public IEnumerable<Point> GenerateGrid(double maxX, double maxY, double zoom, int r)
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
                    if (!(((i) * (i) + (j+ maxY / 4) * (j + maxY / 4) < r * r) || ((i) * (i) + (j - maxY / 4) * (j - maxY / 4) < r * r)))
                    {
                        var point = new Point(i * zoom, j * zoom);

                        if (i == (maxX / 2) || i == -(maxX / 2) || j == -(maxY / 2) || j == (maxY / 2))
                            point.IsBoundary = true;
                        else if ((((i) * (i) + (j + maxY / 4) * (j + maxY / 4) <= (r + 0.25) * (r + 0.25))))
                        {
                            point.IsInnerBoundaryPlus = true;
                            point.IsBoundary = true;
                        }
                        else if (((i) * (i) + (j - maxY / 4) * (j - maxY / 4) <= (r + 0.25) * (r + 0.25)))
                        {
                            point.IsInnerBoundaryMinus = true;
                            point.IsBoundary = true;

                        }
                        points.Add(point);
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


        public void GetAB(List<Point> points, out List<List<double>> aList, out List<double> bList, double plus, double minus)
        {
            aList = new List<List<double>>();
            bList = new List<double>();
            var i = 0;
            foreach (var iPoint in points)
            {
                var b = double.Epsilon;
                //var b =0.0d;
                aList.Add(new List<double>());
                foreach (var jPoint in points)
                {
                    var a = double.Epsilon;
                    //var a = 0.0d;
                    if (!iPoint.IsBoundary && jPoint.IsBoundary)
                    {
                        var commonTriangles = GetCommonTriangles(iPoint, jPoint).ToList();
                        //var borderTriangle = commonTriangles.Aggregate((i1, i2) => i1.BoundaryPointsCount > i2.BoundaryPointsCount ? i1 : i2);
                        //var innerTriangle = commonTriangles.Aggregate((i1, i2) => i1.BoundaryPointsCount < i2.BoundaryPointsCount ? i1 : i2);
                        if (commonTriangles.Count > 1)
                        {
                            var lowerVertex1 = commonTriangles[0].Vertices.Where(coord => coord != iPoint && coord != jPoint).ToList().LastOrDefault();
                            var lowerVertex2 = commonTriangles[1].Vertices.Where(coord => coord != iPoint && coord != jPoint).ToList().LastOrDefault();
                            var lowV1 = new Point3D(lowerVertex1.X, lowerVertex1.Y, 0);
                            var lowV2 = new Point3D(lowerVertex2.X, lowerVertex2.Y, 0);
                            var in0 = new Point3D(jPoint.X, jPoint.Y, 0);
                            var in1 = new Point3D(jPoint.X, jPoint.Y, 1);
                            var bor1 = new Point3D(iPoint.X, iPoint.Y, 1);
                            var bor0 = new Point3D(iPoint.X, iPoint.Y, 0);
                            Triangle3D t11 = new Triangle3D(bor1, lowV1, in0);
                            Triangle3D t12 = new Triangle3D(bor0, lowV1, in1);
                            Triangle3D t21 = new Triangle3D(bor1, lowV2, in0);
                            Triangle3D t22 = new Triangle3D(bor0, lowV2, in1);
                            double c;
                            if (jPoint.IsInnerBoundaryPlus)
                            {
                                c = plus;
                            }
                            else if (jPoint.IsInnerBoundaryMinus)
                            {
                                c = minus;
                            }
                            else
                            {
                                c = double.Epsilon;
                                
                            }
                            
                            
                            b += c*((t11.A * t12.A + t11.B * t12.B) * t11.S0 + (t21.A * t22.A + t21.B * t22.B) * t21.S0);
                        }
                    }
                    else if (!iPoint.IsBoundary && iPoint == jPoint)
                    {
                        var triangles = iPoint.AdjacentTriangles.ToList();
                        foreach (var triangle in triangles)
                        {
                            var lowerVertex = triangle.Vertices.Where(coord => coord != iPoint).ToArray();
                            var m1 = new Point3D(lowerVertex[0].X, lowerVertex[0].Y, 0);
                            var m2 = new Point3D(lowerVertex[1].X, lowerVertex[1].Y, 0);
                            var m3 = new Point3D(iPoint.X, iPoint.Y, 1);
                            Triangle3D t3D = new Triangle3D(m1, m2, m3);
                            a -= (t3D.A * t3D.A + t3D.B * t3D.B) * t3D.S0;
                        }
                    }
                    else if (GetCommonTriangles(iPoint, jPoint).ToList().Count > 0 && !iPoint.IsBoundary && !jPoint.IsBoundary)
                    {
                        var commonTriangles = GetCommonTriangles(iPoint, jPoint).ToList();
                        if (commonTriangles.Count > 1)
                        {
                            var lowerVertex1 = commonTriangles[0].Vertices.Where(coord => coord != iPoint && coord != jPoint).ToList().LastOrDefault();
                            var lowerVertex2 = commonTriangles[1].Vertices.Where(coord => coord != iPoint && coord != jPoint).ToList().LastOrDefault();
                            var lowV1 = new Point3D(lowerVertex1.X, lowerVertex1.Y, 0);
                            var lowV2 = new Point3D(lowerVertex2.X, lowerVertex2.Y, 0);
                            var in0 = new Point3D(jPoint.X, jPoint.Y, 0);
                            var in1 = new Point3D(jPoint.X, jPoint.Y, 1);
                            var bor1 = new Point3D(iPoint.X, iPoint.Y, 1);
                            var bor0 = new Point3D(iPoint.X, iPoint.Y, 0);
                            Triangle3D t11 = new Triangle3D(bor1, lowV1, in0);
                            Triangle3D t12 = new Triangle3D(bor0, lowV1, in1);
                            Triangle3D t21 = new Triangle3D(bor1, lowV2, in0);
                            Triangle3D t22 = new Triangle3D(bor0, lowV2, in1);
                            a -= (t11.A * t12.A + t11.B * t12.B) * t11.S0 + (t21.A * t22.A + t21.B * t22.B) * t21.S0;
                        }
                    }

                    aList[i].Add(a);
                    
                }
                ++i;
                bList.Add(b);
            }
        }
        
        public List<double> DoKachmarz(List<List<double>> a, List<double> b)
        {
            // nn - количество неизвестных;  ny - количество уравнений
            double eps = 1.0e-6;
            //float s;
            int i, j, k;
            double s1, s2, fa1, t;
            double[] x1 = new double[b.Count];
            List<double> x = new List<double>();
            x.Add(double.Epsilon);
            for (i = 1; i < b.Count; i++) x.Add(double.Epsilon);
           

            s1 = s2 = 1d;
            while (s1 > eps * s2)
            {
                for (i = 0; i < b.Count; i++) x1[i] = x[i];

                for (i = 0; i < b.Count; i++)
                {
                    s1 = double.Epsilon;
                    s2 = double.Epsilon;

                    for (j = 0; j < b.Count; j++)
                    {
                        fa1 = a[i][j];
                        s1 += fa1 * x[j];
                        s2 += fa1 * fa1;
                    }
                    t = Convert.ToDouble((double)(b[i] - s1) / s2);
                    if (double.IsInfinity(t))
                        t = (double)1000000;
                    for (k = 0; k < b.Count; k++)
                        x[k] += a[i][k] * t;
                }

                s1 = double.Epsilon;
                s2 = double.Epsilon;

                for (i = 0; i < b.Count; i++)
                {
                    s1 += (x[i] - x1[i]) * (x[i] - x1[i]);
                    s2 += x[i] * x[i];
                }
                s1 = Math.Sqrt(s1);
                s2 = Math.Sqrt(s2);
            }

            return x;
        }

        public List<List<Point>> GetIsolines(int n,  List<Triangle> triangles)
        {
            var z = triangles.Select(item => item.Vertices).ToList();
            var z1 = z.SelectMany(item => item).ToList();
            var maxFi = z1.Max(item => item.Fi);
            var minFi = z1.Min(item => item.Fi);
            var dFi = (maxFi - minFi) / (double)(n+1);
            List<List<Point>> isolines = new List<List<Point>>();
            for (int i = 0; i < n; i++)
            {
               isolines.Add(new List<Point>());
                var fiLevel = minFi + (double)(i+1) * dFi;
                foreach (var triangle in triangles)
                {
                    var fiList = triangle.Vertices.ToList();
                    fiList.Sort((p1, p2) => p1.Fi.CompareTo(p2.Fi));
                    if (!(fiLevel - fiList[2].Fi <0.2 ) || !(fiLevel - fiList[0].Fi >0.2)) continue;
                    double x,y;
                    if (fiLevel - fiList[1].Fi < 0.2)
                    {
                        x = fiList[1].X - (fiList[1].Fi - fiLevel) / (fiList[1].Fi - fiList[0].Fi) * (fiList[1].X - fiList[0].X);
                        y = fiList[1].Y - (fiList[1].Fi - fiLevel) / (fiList[1].Fi - fiList[0].Fi) * (fiList[1].Y - fiList[0].Y);
                    }
                    else
                    {
                        x = fiList[2].X - (fiList[2].Fi - fiLevel) / (fiList[2].Fi - fiList[1].Fi) * (fiList[2].X - fiList[1].X);
                        y = fiList[2].Y - (fiList[2].Fi - fiLevel) / (fiList[2].Fi - fiList[1].Fi) * (fiList[2].Y - fiList[1].Y);
                    }
                    isolines[i].Add(new Point(x,y,fiLevel));

                }

                var noduplicates = isolines[i].Distinct(new PointComparer()).ToList();
                var newIsoline = new List<Point>();
                newIsoline.Add(noduplicates.First());
                foreach (var iPoint in noduplicates)
                {
                    Point point = new Point(0,0);
                    double r = double.MaxValue;
                    foreach (var jPoint in noduplicates)
                    {
                        
                        if (iPoint != jPoint)
                        {

                            if (Math.Sqrt((iPoint.X - jPoint.X) * (iPoint.X - jPoint.X) + (iPoint.Y - jPoint.Y) * (iPoint.Y - jPoint.Y)) < r && !newIsoline.Contains(jPoint))
                            {
                                r = Math.Sqrt((iPoint.X - jPoint.X) * (iPoint.X - jPoint.X) + (iPoint.Y - jPoint.Y) * (iPoint.Y - jPoint.Y));
                                point = jPoint;
                            }
                        }
                    }
                    if(newIsoline.Count<noduplicates.Count && point.X!=0&&point.Y!=0)
                        newIsoline.Add(point);
                }
                isolines[i] = newIsoline;
            }

            return isolines;
        }

        public List<List<Point>> GetFieldlines(List<Triangle> triangles, List<Point> points, double maxX, double maxY, double fiMax, double coef)
        {
            List<List<Point>> fieldlines = new List<List<Point>>();
            var boundaryPoints = triangles.Where(tri => tri.SignBoundaryPointsCount>0).Select(tri => tri.Circumcenter).ToList();
            var usedTriangles = new List<Triangle>();
            int i = 0;
            foreach (var point in boundaryPoints)
            {
                Triangle t2d = new Triangle(new Point(0,10), new Point(10,20), point );
                Triangle3D t3d = new Triangle3D(new Point3D(), new Point3D(), new Point3D());
                var p1 = new Point3D(t2d.Vertices[0].X, t2d.Vertices[0].Y, t2d.Vertices[0].Fi);
                var p2 = new Point3D(t2d.Vertices[1].X, t2d.Vertices[1].Y, t2d.Vertices[1].Fi);
                var p3 = new Point3D(t2d.Vertices[2].X, t2d.Vertices[2].Y, t2d.Vertices[2].Fi);
                fieldlines.Add(new List<Point>());
                foreach (var triangle in triangles)
                {
                    if (!point.PointInTriangle(point, triangle.Vertices[0], triangle.Vertices[1], triangle.Vertices[2])) continue;
                    t2d = triangle;
                    usedTriangles.Add(t2d);
                    p1 = new Point3D(t2d.Vertices[0].X, t2d.Vertices[0].Y, t2d.Vertices[0].Fi);
                    p2 = new Point3D(t2d.Vertices[1].X, t2d.Vertices[1].Y, t2d.Vertices[1].Fi);
                    p3 = new Point3D(t2d.Vertices[2].X, t2d.Vertices[2].Y, t2d.Vertices[2].Fi);

                    t3d = new Triangle3D(p1, p2, p3);
                    break;
                }

                bool intr = true;
                var curX = point.X;
                var curY = point.Y;
                var curFi = -(t3d.A*curX + t3d.B * curY + t3d.D)/ t3d.C;
                while (curX<=maxX/2&&curX>=-maxX/2&& curY <= maxY / 2 && curY >= -maxY / 2 && Math.Abs(curFi) <=fiMax*1.5)
                {
                    var addPoint = new Point(curX, curY, curFi);
                    fieldlines[i].Add(addPoint);
                    if(!intr)
                        break;
                    if (addPoint.PointInTriangle(addPoint, t2d.Vertices[0], t2d.Vertices[1], t2d.Vertices[2]))
                    {
                        Vector3D v1 = p2 - p1;
                        Vector3D v2 = p3 - p2;

                        // Get the cross product.
                        Vector3D n = Vector3D.CrossProduct(v1, v2);
                        n.Normalize();
                        if (curFi < 0)
                        {
                            n.X = -n.X;
                            n.Y = -n.Y;
                            n.Z = -n.Z;
                        }

                        curX += coef * n.X;
                        curY += coef * n.Y;
                        curFi += coef * n.Z;

                        //if (t3d.A * t3d.D < 0)
                        //    curX += coef * Math.Abs(t3d.A);
                        //else
                        //    curX -= coef * Math.Abs(t3d.A);
                        //if (t3d.B * t3d.D < 0)
                        //    curY += coef * Math.Abs(t3d.B);
                        //else
                        //    curY -= coef * Math.Abs(t3d.B);
                        //if (t3d.D > 0)
                        //    curFi -= coef * Math.Abs(t3d.C);
                        //else
                        //    curFi += coef * Math.Abs(t3d.C);

                        //if(coef * Math.Abs(t3d.C) < 0.00001)
                        //    break;
                    }
                    else
                    {
                        intr = false;
                        foreach (var triangle in triangles)
                        {
                            if (addPoint.PointInTriangle(addPoint, triangle.Vertices[0], triangle.Vertices[1], triangle.Vertices[2]))
                            {
                                t2d = triangle;
                                p1 = new Point3D(t2d.Vertices[0].X, t2d.Vertices[0].Y, t2d.Vertices[0].Fi);
                                p2 = new Point3D(t2d.Vertices[1].X, t2d.Vertices[1].Y, t2d.Vertices[1].Fi);
                                p3 = new Point3D(t2d.Vertices[2].X, t2d.Vertices[2].Y, t2d.Vertices[2].Fi);
                                t3d = new Triangle3D(p1, p2, p3);
                                intr = true;
                                break;
                            }
                        }
                    }


                }

                ++i;
            }

            return fieldlines;
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
