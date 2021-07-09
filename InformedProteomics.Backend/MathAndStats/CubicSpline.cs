//
// Author: Ryan Seghers
//
// Copyright (C) 2013-2014 Ryan Seghers
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the irrevocable, perpetual, worldwide, and royalty-free
// rights to use, copy, modify, merge, publish, distribute, sublicense,
// display, perform, create derivative works from and/or sell copies of
// the Software, both in source and object code form, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

using System;
using PRISM;

namespace InformedProteomics.Backend.MathAndStats
{
    /// <summary>
    /// Cubic spline interpolation
    /// <para>
    /// Call Fit (or use the corrector constructor) to compute spline coefficients, then Eval to evaluate the spline at other X coordinates
    /// </para>
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is implemented based on the Wikipedia article:
    /// http://en.wikipedia.org/wiki/Spline_interpolation
    /// I'm not sure I have the right to include a copy of the article so the equation numbers referenced in
    /// comments will end up being wrong at some point.
    /// </para>
    /// <para>
    /// This is not optimized, and is not MT safe.
    /// This can extrapolate off the ends of the splines.
    /// You must provide points in X sort order.
    /// </para>
    /// </remarks>
    public class CubicSpline
    {
        // Ignore Spelling: eval, tri, Wikipedia, xs

        #region Fields

        // N-1 spline coefficients for N points
        private float[] a;
        private float[] b;

        // Save the original x and y for Eval
        private float[] xOrig;
        private float[] yOrig;

        #endregion

        #region Constructor

        /// <summary>
        /// Constructor
        /// </summary>
        public CubicSpline()
        {
        }

        /// <summary>
        /// Construct and call Fit
        /// </summary>
        /// <param name="x">Input: X coordinates to fit</param>
        /// <param name="y">Input: Y coordinates to fit</param>
        /// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint</param>
        /// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint</param>
        /// <param name="debug">True to write messages to the console</param>
        public CubicSpline(float[] x, float[] y, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
        {
            Fit(x, y, startSlope, endSlope, debug);
        }

        #endregion

        #region Private Methods

        /// <summary>
        /// Throws an exception if Fit has not yet been called
        /// </summary>
        private void CheckAlreadyFitted()
        {
            if (a == null)
            {
                throw new Exception("Fit must be called before you can evaluate.");
            }
        }

        private int _lastIndex;

        /// <summary>
        /// Find where in xOrig the specified x falls, by simultaneous traverse.
        /// This allows xs to be less than x[0] and/or greater than x[n-1]. So allows extrapolation.
        /// This keeps state, so requires that x be sorted and xs called in ascending order, and is not multi-thread safe.
        /// </summary>
        private int GetNextXIndex(float x)
        {
            if (x < xOrig[_lastIndex])
            {
                throw new ArgumentException("The X values to evaluate must be sorted.");
            }

            while ((_lastIndex < xOrig.Length - 2) && (x > xOrig[_lastIndex + 1]))
            {
                _lastIndex++;
            }

            return _lastIndex;
        }

        /// <summary>
        /// Evaluate the specified x value using the specified spline
        /// </summary>
        /// <param name="x">The x value</param>
        /// <param name="j">Which spline to use</param>
        /// <param name="debug">True to enable writing messages to the console</param>
        /// <returns>The y value</returns>
        private float EvalSpline(float x, int j, bool debug = false)
        {
            var dx = xOrig[j + 1] - xOrig[j];
            var t = (x - xOrig[j]) / dx;
            var y = (1 - t) * yOrig[j] + t * yOrig[j + 1] + t * (1 - t) * (a[j] * (1 - t) + b[j] * t); // equation 9
            if (debug)
            {
                ShowDebug(string.Format("xs = {0}, j = {1}, t = {2}", x, j, t));
            }

            return y;
        }

        #endregion

        #region Fit*

        /// <summary>
        /// Fit x,y and then evaluate at points in xs and return the corresponding y's.
        /// This does the "natural spline" style for ends.
        /// This can extrapolate off the ends of the splines.
        /// You must provide points in X sort order.
        /// </summary>
        /// <param name="x">Input: X coordinates to fit</param>
        /// <param name="y">Input: Y coordinates to fit</param>
        /// <param name="xs">Input: X coordinates to evaluate the fitted curve at</param>
        /// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint</param>
        /// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint</param>
        /// <param name="debug">Turn on console output. Default is false</param>
        /// <returns>The computed y values for each xs</returns>
        public float[] FitAndEval(float[] x, float[] y, float[] xs, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
        {
            Fit(x, y, startSlope, endSlope, debug);
            return Eval(xs, debug);
        }

        /// <summary>
        /// Compute spline coefficients for the specified x,y points.
        /// This does the "natural spline" style for ends.
        /// This can extrapolate off the ends of the splines.
        /// You must provide points in X sort order.
        /// </summary>
        /// <param name="x">Input: X coordinates to fit</param>
        /// <param name="y">Input: Y coordinates to fit</param>
        /// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint</param>
        /// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint</param>
        /// <param name="debug">Turn on console output. Default is false</param>
        public void Fit(float[] x, float[] y, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
        {
            if (float.IsInfinity(startSlope) || float.IsInfinity(endSlope))
            {
                throw new Exception("startSlope and endSlope cannot be infinity.");
            }

            // Save x and y for eval
            xOrig = x;
            yOrig = y;

            var n = x.Length;
            var r = new float[n]; // the right hand side numbers: Wikipedia page overloads b

            var m = new TriDiagonalMatrixF(n);
            float dx1;
            float dy1;

            // First row is different (equation 16 from the article)
            if (float.IsNaN(startSlope))
            {
                dx1 = x[1] - x[0];
                m.C[0] = 1.0f / dx1;
                m.B[0] = 2.0f * m.C[0];
                r[0] = 3 * (y[1] - y[0]) / (dx1 * dx1);
            }
            else
            {
                m.B[0] = 1;
                r[0] = startSlope;
            }

            // Body rows (equation 15 from the article)
            for (var i = 1; i < n - 1; i++)
            {
                dx1 = x[i] - x[i - 1];
                var dx2 = x[i + 1] - x[i];

                m.A[i] = 1.0f / dx1;
                m.C[i] = 1.0f / dx2;
                m.B[i] = 2.0f * (m.A[i] + m.C[i]);

                dy1 = y[i] - y[i - 1];
                var dy2 = y[i + 1] - y[i];
                r[i] = 3 * (dy1 / (dx1 * dx1) + dy2 / (dx2 * dx2));
            }

            // Last row also different (equation 17 from the article)
            if (float.IsNaN(endSlope))
            {
                dx1 = x[n - 1] - x[n - 2];
                dy1 = y[n - 1] - y[n - 2];
                m.A[n - 1] = 1.0f / dx1;
                m.B[n - 1] = 2.0f * m.A[n - 1];
                r[n - 1] = 3 * (dy1 / (dx1 * dx1));
            }
            else
            {
                m.B[n - 1] = 1;
                r[n - 1] = endSlope;
            }

            if (debug)
            {
                ShowDebug(string.Format("Tri-diagonal matrix:\n{0}", m.ToDisplayString(":0.0000", "  ")));
            }

            if (debug)
            {
                ShowDebug(string.Format("r: {0}", ArrayUtil.ToString(r)));
            }

            // k is the solution to the matrix
            var k = m.Solve(r);
            if (debug)
            {
                ShowDebug(string.Format("k = {0}", ArrayUtil.ToString(k)));
            }

            // a and b are each spline's coefficients
            a = new float[n - 1];
            b = new float[n - 1];

            for (var i = 1; i < n; i++)
            {
                dx1 = x[i] - x[i - 1];
                dy1 = y[i] - y[i - 1];
                a[i - 1] = k[i - 1] * dx1 - dy1; // equation 10 from the article
                b[i - 1] = -k[i] * dx1 + dy1; // equation 11 from the article
            }

            if (debug)
            {
                ShowDebug(string.Format("a: {0}", ArrayUtil.ToString(a)));
            }

            if (debug)
            {
                ShowDebug(string.Format("b: {0}", ArrayUtil.ToString(b)));
            }
        }

        #endregion

        #region Eval*

        /// <summary>
        /// Evaluate the spline at the specified x coordinates.
        /// This can extrapolate off the ends of the splines.
        /// You must provide X's in ascending order.
        /// The spline must already be computed before calling this, meaning you must have already called Fit() or FitAndEval().
        /// </summary>
        /// <param name="x">Input. X coordinates to evaluate the fitted curve at</param>
        /// <param name="debug">Turn on console output. Default is false</param>
        /// <returns>The computed y values for each x</returns>
        public float[] Eval(float[] x, bool debug = false)
        {
            CheckAlreadyFitted();

            var n = x.Length;
            var y = new float[n];
            _lastIndex = 0; // Reset simultaneous traversal in case there are multiple calls

            for (var i = 0; i < n; i++)
            {
                // Find which spline can be used to compute this x (by simultaneous traverse)
                var j = GetNextXIndex(x[i]);

                // Evaluate using j'th spline
                y[i] = EvalSpline(x[i], j, debug);
            }

            return y;
        }

        /// <summary>
        /// Evaluate (compute) the slope of the spline at the specified x coordinates.
        /// This can extrapolate off the ends of the splines.
        /// You must provide X's in ascending order.
        /// The spline must already be computed before calling this, meaning you must have already called Fit() or FitAndEval().
        /// </summary>
        /// <param name="x">Input. X coordinates to evaluate the fitted curve at</param>
        /// <param name="debug">Turn on console output. Default is false</param>
        /// <returns>The computed y values for each x</returns>
        public float[] EvalSlope(float[] x, bool debug = false)
        {
            CheckAlreadyFitted();

            var n = x.Length;
            var qPrime = new float[n];
            _lastIndex = 0; // Reset simultaneous traversal in case there are multiple calls

            for (var i = 0; i < n; i++)
            {
                // Find which spline can be used to compute this x (by simultaneous traverse)
                var j = GetNextXIndex(x[i]);

                // Evaluate using j'th spline
                var dx = xOrig[j + 1] - xOrig[j];
                var dy = yOrig[j + 1] - yOrig[j];
                var t = (x[i] - xOrig[j]) / dx;

                // From equation 5 we could also compute q' (qp) which is the slope at this x
                qPrime[i] = dy / dx
                    + (1 - 2 * t) * (a[j] * (1 - t) + b[j] * t) / dx
                    + t * (1 - t) * (b[j] - a[j]) / dx;

                if (debug)
                {
                    ShowDebug(string.Format("[{0}]: xs = {1}, j = {2}, t = {3}", i, x[i], j, t));
                }
            }

            return qPrime;
        }

        private void ShowDebug(string message)
        {
            ConsoleMsgUtils.ShowDebug(message);
        }

        #endregion

        #region Static Methods

        /// <summary>
        /// Static all-in-one method to fit the splines and evaluate at X coordinates
        /// </summary>
        /// <param name="x">Input: X coordinates to fit</param>
        /// <param name="y">Input: Y coordinates to fit</param>
        /// <param name="xs">Input: X coordinates to evaluate the fitted curve at</param>
        /// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint</param>
        /// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint</param>
        /// <param name="debug">Turn on console output. Default is false</param>
        /// <returns>The computed y values for each xs</returns>
        public static float[] Compute(float[] x, float[] y, float[] xs, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
        {
            var spline = new CubicSpline();
            return spline.FitAndEval(x, y, xs, startSlope, endSlope, debug);
        }

        /// <summary>
        /// Fit the input x,y points using a 'geometric' strategy so that y does not have to be a single-valued
        /// function of x
        /// </summary>
        /// <param name="x">Input: x coordinates</param>
        /// <param name="y">Input: y coordinates, do not need to be a single-valued function of x</param>
        /// <param name="nOutputPoints">How many output points to create</param>
        /// <param name="xs">Output (interpolated) x values</param>
        /// <param name="ys">Output (interpolated) y values</param>
        public static void FitGeometric(float[] x, float[] y, int nOutputPoints, out float[] xs, out float[] ys)
        {
            // Compute distances
            var n = x.Length;
            var cumulativeDistances = new float[n];
            cumulativeDistances[0] = 0;
            float totalDist = 0;

            for (var i = 1; i < n; i++)
            {
                var dx = x[i] - x[i - 1];
                var dy = y[i] - y[i - 1];
                var dist = (float)Math.Sqrt(dx * dx + dy * dy);
                totalDist += dist;
                cumulativeDistances[i] = totalDist;
            }

            // Create 'times' to interpolate to
            var dt = totalDist / (nOutputPoints - 1);
            var times = new float[nOutputPoints];
            times[0] = 0;

            for (var i = 1; i < nOutputPoints; i++)
            {
                times[i] = times[i - 1] + dt;
            }

            // Spline fit both x and y to times
            var xSpline = new CubicSpline();
            xs = xSpline.FitAndEval(cumulativeDistances, x, times);

            var ySpline = new CubicSpline();
            ys = ySpline.FitAndEval(cumulativeDistances, y, times);
        }

        #endregion
    }
}
