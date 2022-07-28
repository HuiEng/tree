// ************** Generate the tree diagram	 *****************
var margin = { top: 20, right: 120, bottom: 20, left: 120 },
	width = 3000 - margin.right - margin.left,
	height = 2000 - margin.top - margin.bottom,
	branch_len = 50;

var i = 0,
	duration = 750,
	root;

var tree = d3.layout.tree()
	.size([height, width]);

var diagonal = d3.svg.diagonal()
	.projection(function (d) { return [d.x, d.y]; });

var svg = d3.select("body").append("svg")
	.attr("width", width + margin.right + margin.left)
	.attr("height", height + margin.top + margin.bottom)
	.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

root = treeData;
root.x0 = height / 2;
root.y0 = 0;

update(root);

d3.select(self.frameElement).style("height", "500px");

function wrap(text, width) {
	text.each(function () {
		var text = d3.select(this),
			words = text.text().split(/,/).reverse(),
			word,
			line = [],
			lineNumber = 0,
			lineHeight = 0.02, // ems
			x = text.attr("x"),
			y = text.attr("y"),
			dy = parseFloat(text.attr("dy")),
			tspan = text.text(null)
				.append("tspan")
				.attr("x", x)
				.attr("y", y)
				.attr("dy", dy + "em");
		while (word = words.pop()) {
			line.push(word);
			tspan.text(line.join(" "));
			if (tspan.node().getComputedTextLength() > width) {
				line.pop();
				tspan.text(line.join(" "));
				line = [word];
				tspan = text.append("tspan")
					.attr("x", x)
					.attr("y", y)
					.attr("dy", ++lineNumber * lineHeight + dy + "em")
					.text(word);
			}
		}
	});
}



function update(source) {

	// Compute the new tree layout.
	var nodes = tree.nodes(root).reverse(),
		links = tree.links(nodes);

	// Normalize for fixed-depth.
	nodes.forEach(function (d) { d.y = d.depth * branch_len; });

	// Update the nodes…
	var node = svg.selectAll("g.node")
		.data(nodes, function (d) { return d.id || (d.id = ++i); });

	// Enter any new nodes at the parent's previous position.
	var nodeEnter = node.enter().append("g")
		.attr("class", "node")
		.attr("transform", function (d) { return "translate(" + source.x0 + "," + source.y0 + ")"; })
		.attr("children", function (d) {
			if (d.children) {
				// if (d.children.length > 5) {
				// 	d._children = d.children;
				// 	d.children = null;
				// }
				// if (parseInt(d.node) <1230) {
				// 	d._children = d.children;
				// 	d.children = null;
				// }
				// console.log(d.node);
				// var child = d.children;
				// while (child.children){
				// 	child = child.children;
				// 	console.log(child);
				// }
				// if (child[0].node == ">>>") {
				// 	d._children = d.children;
				// 	d.children = null;
				// }
			}
		})
		.on("click", click);
	//   .on("mouseover",mouseover);

	nodeEnter.append("circle")
		.attr("r", 1e-6)
		//.style("fill", function (d) { {return d.node ? "lightsteelblue" : "#fff"; }});
		// .style("fill", function (d) { return d._children ? "lightsteelblue" : "#fff"; });
		.style("fill", function (d) {return d.content.split(/,/).length == 2 ? "red" : "lightsteelblue";});

	nodeEnter.append("text")
		.attr("x", function (d) { return d.children || d._children ? -13 : 1; })
		.attr("dy", "1.5em")
		.attr("text-anchor", "middle")
		.text(function (d) {
			// if (d.children) {
			// 	return "";
			// } else {
			// 	return d.node;
			// }
			//?
			//return d.node - 1;
			// return d.childCount;
			return d.content;
		}).call(wrap, 1)
		.style("fill-opacity", 1e-6);

	// Transition nodes to their new position.
	var nodeUpdate = node.transition()
		.duration(duration)
		.attr("transform", function (d) { return "translate(" + d.x + "," + d.y + ")"; });

	nodeUpdate.select("circle")
		.attr("r", 10)
		// .style("fill", function (d) { return d._children ? "lightsteelblue" : "#fff"; })
		.style("fill", function (d) {return d.content.split(/,/).length == 2 ? "red" : "lightsteelblue";})
		.style("stroke", function (d) {return d.content.split(/,/).length == 2 ? "red" : "lightsteelblue";});
		// .style("stroke", function (d) { return d.node ? "lightsteelblue" : "#fff"; });

	nodeUpdate.select("text")
		.style("fill-opacity", 1);

	// Transition exiting nodes to the parent's new position.
	var nodeExit = node.exit().transition()
		.duration(duration)
		.attr("transform", function (d) { return "translate(" + source.x + "," + source.y + ")"; })
		.remove();

	nodeExit.select("circle")
		.attr("r", 1e-6);

	nodeExit.select("text")
		.style("fill-opacity", 1e-6);

	// Update the links…
	var link = svg.selectAll("path.link")
		.data(links, function (d) { if (d.target.node) { return d.target.id; } });

	// Enter any new links at the parent's previous position.
	link.enter().insert("path", "g")
		.attr("class", "link")
		.attr("d", function (d) {
			var o = { x: source.x0, y: source.y0 };
			return diagonal({ source: o, target: o });
		});

	// Transition links to their new position.
	link.transition()
		.duration(duration)
		.attr("d", diagonal);

	// Transition exiting nodes to the parent's new position.
	link.exit().transition()
		.duration(duration)
		.attr("d", function (d) {
			var o = { x: source.x, y: source.y };
			return diagonal({ source: o, target: o });
		})
		.remove();

	// Stash the old positions for transition.
	nodes.forEach(function (d) {
		d.x0 = d.x;
		d.y0 = d.y;
	});
}

// Toggle children on click.
function click(d) {
	//   if (d.children) {
	// 	d._children = d.children;
	// 	d.children = null;
	//   } else {
	// 	d.children = d._children;
	// 	d._children = null;
	//   }
	// update(d);
	mouseover(d);
}

function mouseover(d) {
	document.getElementById("info").innerHTML = "node: " + d.node + "\ncontent: " + d.content + "\n";
}

// click(root)