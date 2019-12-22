#!/usr/bin/env python3

from argparse import ArgumentParser, RawTextHelpFormatter
import collections
import fnmatch
import folium
from folium.plugins import HeatMap
import ijson
import json
import os
from progressbar import ProgressBar, Bar, ETA, Percentage
from utils import *
import webbrowser
from xml.etree import ElementTree
from xml.dom import minidom
import zipfile
from pygeodesy.sphericalTrigonometry import LatLon


class Generator:
    def __init__(self):
        self.coordinates = collections.defaultdict(int)
        self.max_coordinates = (0, 0)
        self.max_magnitude = 0

        self.flight_coordinates = []

    def loadJSONData(self, json_file, date_range):
        """Loads the Google location data from the given json file.

        Arguments:
            json_file {file} -- An open file-like object with JSON-encoded
                Google location data.
            date_range {tuple} -- A tuple containing the min-date and max-date.
                e.g.: (None, None), (None, '2019-01-01'), ('2017-02-11'), ('2019-01-01')
        """
        data = json.load(json_file)
        w = [Bar(), Percentage(), " ", ETA()]
        with ProgressBar(max_value=len(data["locations"]), widgets=w) as pb:
            for i, loc in enumerate(data["locations"]):
                if "latitudeE7" not in loc or "longitudeE7" not in loc:
                    continue
                coords = (round(loc["latitudeE7"] / 1e7, 6),
                           round(loc["longitudeE7"] / 1e7, 6))

                if timestampInRange(loc['timestampMs'], date_range):
                    self.updateCoord(coords)
                pb.update(i)

    def streamJSONData(self, json_file, date_range):
        """Stream the Google location data from the given json file.

        Arguments:
            json_file {file} -- An open file-like object with JSON-encoded
                Google location data.
            date_range {tuple} -- A tuple containing the min-date and max-date.
                e.g.: (None, None), (None, '2019-01-01'), ('2017-02-11'), ('2019-01-01')
        """
        # Estimate location amount
        max_value_est = sum(1 for line in json_file) / 13
        json_file.seek(0)

        locations = ijson.items(json_file, "locations.item")
        w = [Bar(), Percentage(), " ", ETA()]
        with ProgressBar(max_value=max_value_est, widgets=w) as pb:
            for i, loc in enumerate(locations):
                if "latitudeE7" not in loc or "longitudeE7" not in loc:
                    continue
                coords = (round(loc["latitudeE7"] / 1e7, 6),
                            round(loc["longitudeE7"] / 1e7, 6))

                if timestampInRange(loc['timestampMs'], date_range):
                    self.updateCoord(coords)

                if i > max_value_est:
                    max_value_est = i
                    pb.max_value = i
                pb.update(i)

    def loadKMLData(self, file_name, date_range):
        """Loads the Google location data from the given KML file.

        Arguments:
            file_name {string or file} -- The name of the KML file
                (or an open file-like object) with the Google location data.
            date_range {tuple} -- A tuple containing the min-date and max-date.
                e.g.: (None, None), (None, '2019-01-01'), ('2017-02-11'), ('2019-01-01')
        """
        xmldoc = minidom.parse(file_name)
        gxtrack = xmldoc.getElementsByTagName("gx:coord")
        when = xmldoc.getElementsByTagName("when")
        w = [Bar(), Percentage(), " ", ETA()]

        with ProgressBar(max_value=len(gxtrack), widgets=w) as pb:
            for i, number in enumerate(gxtrack):
                loc = (number.firstChild.data).split()
                coords = (round(float(loc[1]), 6), round(float(loc[0]), 6))
                date = when[i].firstChild.data
                if dateInRange(date[:10], date_range):
                    self.updateCoord(coords)
                pb.update(i)

    def loadFlightsData(self, flights_path, date_range):
        """Loads the Google semantic location data from the given folder.

        Arguments:
            flights_path {string} -- Semantic Location History folder from
                Google Takeout Zip archive.
            date_range {tuple} -- A tuple containing the min-date and max-date.
                e.g.: (None, None), (None, '2019-01-01'), ('2017-02-11'), ('2019-01-01')
        """
        for year_folder in os.listdir(flights_path):
            year_path = os.path.join(flights_path, year_folder)
            if not os.path.isdir(year_path):
                continue

            for month_file in os.listdir(year_path):
                if not month_file.endswith(".json"):
                    continue

                month_path = os.path.join(year_path, month_file)
                with open(month_path) as json_file:
                    data = json.load(json_file)
                    if "timelineObjects" not in data:
                        continue

                    for i, loc in enumerate(data["timelineObjects"]):
                        if "activitySegment" not in loc:
                            continue

                        segment = loc["activitySegment"]
                        if "activityType" not in segment or \
                                "confidence" not in segment:
                            continue

                        start_timestamp = segment["duration"]["startTimestampMs"]
                        if segment["activityType"] != "FLYING" or \
                                not (segment["confidence"] == "HIGH" or \
                                     segment["confidence"] == "MEDIUM") or \
                                not timestampInRange(start_timestamp, date_range):
                            continue

                        start = segment["startLocation"]
                        end = segment["endLocation"]

                        start = LatLon(start["latitudeE7"] / 1e7, start["longitudeE7"] / 1e7)
                        end = LatLon(end["latitudeE7"] / 1e7, end["longitudeE7"] / 1e7)

                        intermediates = []
                        for i in range(0, 101):
                            location = start.intermediateTo(end, i / 100.0)
                            intermediates.append([location.lat, location.lon ])

                        self.flight_coordinates.append(intermediates)

    def loadZIPData(self, file_name, date_range):
        # TODO: need to make this work for the flights setting
        """
        Load Google location data from a "takeout-*.zip" file.
        """
        from bs4 import BeautifulSoup
        """
        <div class="service_name">
            <h1 class="data-folder-name" data-english-name="LOCATION_HISTORY" data-folder-name="Location History">
                Location History
            </h1>
        </div>
        """
        zip_file = zipfile.ZipFile(file_name)
        namelist = zip_file.namelist()
        (html_path,) = fnmatch.filter(namelist, "Takeout/*.html")
        with zip_file.open(html_path) as read_file:
            soup = BeautifulSoup(read_file, "html.parser")
        (elem,) = soup.select(
            "#service-tile-LOCATION_HISTORY > button > div.service_summary > div > h1[data-english-name=LOCATION_HISTORY]")
        name = elem["data-folder-name"]
        (data_path,) = fnmatch.filter(
            namelist,
            "Takeout/{name}/{name}.*".format(name=name))
        print("Reading location data file from zip archive: {!r}".format(
            data_path))

        if data_path.endswith(".json"):
            with zip_file.open(data_path) as read_file:
                self.loadJSONData(read_file, date_range)
        elif data_path.endswith(".kml"):
            with zip_file.open(data_path) as read_file:
                self.loadKMLData(read_file, date_range)
        else:
            raise ValueError("unsupported extension for {!r}: only .json and .kml supported"
                .format(file_name))

    def updateCoord(self, coords):
        self.coordinates[coords] += 1
        if self.coordinates[coords] > self.max_magnitude:
            self.max_coordinates = coords
            self.max_magnitude = self.coordinates[coords]

    def generateMap(self, tiles, map_zoom_start=6, heatmap_radius=7,
                    heatmap_blur=4, heatmap_min_opacity=0.2,
                    heatmap_max_zoom=4):
        map_data = [(coords[0], coords[1], magnitude)
                    for coords, magnitude in self.coordinates.items()]

        # Generate map
        map = folium.Map(location=self.max_coordinates,
                       zoom_start=map_zoom_start,
                       tiles=tiles)

        # Generate heat map
        heatmap = HeatMap(map_data,
                          max_val=self.max_magnitude,
                          min_opacity=heatmap_min_opacity,
                          radius=heatmap_radius,
                          blur=heatmap_blur,
                          max_zoom=heatmap_max_zoom)

        map.add_child(heatmap)

        for flight_coordinate in self.flight_coordinates:
            line = folium.PolyLine(flight_coordinate, color="green", weight=2.5, opacity=.25)
            line.add_to(map)

        return map

    def run(self, data_files, output_file, date_range, stream_data, tiles, flights_path):
        """Load the data, generate the heatmap and save it.

        Arguments:
            data_files {list} -- List of names of the data files with the Google
                location data or the Google Takeout ZIP archive.
            output_file {string} -- The name of the output file.
            flights_path {string} -- Semantic Location History folder from
                Google Takeout Zip archive.
        """
        for i, data_file in enumerate(data_files):
            print("\n({}/{}) Loading data from {}".format(
                i + 1,
                len(data_files) + 2,
                data_file))
            if data_file.endswith(".zip"):
                self.loadZIPData(data_file, date_range)
            elif data_file.endswith(".json"):
                with open(data_file) as json_file:
                    if stream_data:
                        self.streamJSONData(json_file, date_range)
                    else:
                        self.loadJSONData(json_file, date_range)
            elif data_file.endswith(".kml"):
                self.loadKMLData(data_file, date_range)
            else:
                raise NotImplementedError(
                    "Unsupported file extension for {!r}".format(data_file))

        total_steps = len(data_files) + 2
        if flights_path is not None:
            total_steps += 1
            print("\n({}/{}) Loading flights data".format(
                total_steps - 2,
                total_steps))
            self.loadFlightsData(flights_path, date_range)

        print("\n({}/{}) Generating heatmap".format(
            total_steps - 1,
            total_steps))
        map = self.generateMap(tiles)
        print("\n({}/{}) Saving map to {}\n".format(
            total_steps,
            total_steps,
            output_file))
        map.save(output_file)


if __name__ == "__main__":
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "files", metavar="file", type=str, nargs='+', help="Any of the following files:\n"
        "1. Your location history JSON file from Google Takeout\n"
        "2. Your location history KML file from Google Takeout\n"
        "3. The takeout-*.zip raw download from Google Takeout \nthat contains either of the above files")
    parser.add_argument("-o", "--output", dest="output", metavar="", type=str, required=False,
                        help="Path of heatmap HTML output file.", default="heatmap.html")
    parser.add_argument("--min-date", dest="min_date", metavar="YYYY-MM-DD", type=str, required=False,
                        help="The earliest date from which you want to see data in the heatmap.")
    parser.add_argument("--max-date", dest="max_date", metavar="YYYY-MM-DD", type=str, required=False,
                        help="The latest date from which you want to see data in the heatmap.")
    parser.add_argument("-s", "--stream", dest="stream", action="store_true", help="Option to iteratively load data.")
    parser.add_argument("--map", "-m", dest="map", metavar="MAP", type=str, required=False, default="OpenStreetMap",
                        help="The name of the map tiles you want to use.\n" \
                        "(e.g. 'OpenStreetMap', 'StamenTerrain', 'StamenToner', 'StamenWatercolor')")
    parser.add_argument("-f", "--flights", dest="flights_path", metavar="", type=str, required=False,
                        help="Path of semantic location history folder " \
                        "(e.g. ./Takeout/Location\ History/Semantic\ Location\ History/)")

    args = parser.parse_args()
    data_file = args.files
    output_file = args.output
    date_range = args.min_date, args.max_date
    tiles = args.map
    stream_data = args.stream
    flights_path = args.flights_path

    generator = Generator()
    generator.run(data_file, output_file, date_range, stream_data, tiles, flights_path)
    # Check if browser is text-based
    if not isTextBasedBrowser(webbrowser.get()):
        try:
            print("[info] Opening {} in browser".format(output_file))
            webbrowser.open("file://" + os.path.realpath(output_file))
        except webbrowser.Error:
            print("[info] No runnable browser found. Open {} manually.".format(
                output_file))
            print("[info] Path to heatmap file: \"{}\"".format(os.path.abspath(output_file)))
