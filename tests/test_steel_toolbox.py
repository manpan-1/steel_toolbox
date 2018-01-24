#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `steel_toolbox` package."""


import unittest

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from click.testing import CliRunner

import steel_toolbox as st
from steel_toolbox import cli


class TestSteel_toolbox(unittest.TestCase):
    """Tests for `steel_toolbox` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_000_something(self):
        """Test something."""

        # The following code is used to cross-check the validity of the results. Points for 2 planes are created and used
        # for testing the fitting and plotting methods.

        # Test xyz data from randomised z=1x+2y+3 for x,y values from -10 t0 +10
        def f_3(beta, xy):
            """ implicit definition of the plane"""
            return beta[0] * xy[0] + beta[1] * xy[1] + beta[2]

        def make_plane_data(beta):
            x = [x + np.random.rand() for x in np.linspace(-10, 10, 21)] + [x + np.random.rand() for x in
                                                                            np.linspace(-10, 10, 21)] + [
                    x + np.random.rand() for x in np.linspace(-10, 10, 21)]
            y = [x + np.random.rand() for x in np.linspace(-10, 10, 21)] + [x + np.random.rand() for x in
                                                                            np.linspace(0, 20, 21)] + [
                    x + np.random.rand()
                    for x in
                    np.linspace(10, 30,
                                21)]
            z = f_3(beta, np.row_stack([x, y]))
            x = np.r_[x]
            y = np.r_[y]
            z = np.r_[z]
            return st.scan_3D.FlatFace(scanned_data=np.transpose(np.row_stack([x, y, z])))

        # Create points for the two planes.
        p1 = make_plane_data([1, 0, -4])
        p2 = make_plane_data([0, 1, -4])

        # Fit a plane on the created points.
        p1.fit_plane()
        p2.fit_plane()
        lp12 = p1.ref_plane & p2.ref_plane

        # Plot the results.
        fig2 = plt.figure()
        Axes3D(fig2)
        p1.plot_face(fig=fig2)
        p2.plot_face(fig=fig2)
        lp12.plot_line(fig=fig2, ends=[-10, 10])

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        assert 'steel_toolbox.cli.main' in result.output
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output
