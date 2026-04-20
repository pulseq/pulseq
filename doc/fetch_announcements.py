#!/usr/bin/env python3
"""Fetch the latest Announcements from GitHub Discussions and generate an
HTML snippet for inclusion in the Doxygen-based project website.

Usage
-----
  # Set a GitHub token (classic PAT with public_repo, or a fine-grained
  # token with read access to Discussions).  In GitHub Actions the
  # built-in GITHUB_TOKEN works.
  export GITHUB_TOKEN="ghp_..."
  python3 fetch_announcements.py            # writes announcements.html
  python3 fetch_announcements.py -n 5       # fetch 5 most recent
  python3 fetch_announcements.py -o out.html

The generated file is meant to be included in a Doxygen .dox page via
  \\htmlinclude announcements.html
"""

import argparse
import json
import os
import sys
import urllib.request
import urllib.error
from datetime import datetime, timezone

GITHUB_GRAPHQL_URL = "https://api.github.com/graphql"
DEFAULT_REPO_OWNER = "pulseq"
DEFAULT_REPO_NAME = "pulseq"
DEFAULT_CATEGORY = "Announcements"
DEFAULT_COUNT = 3


def fetch_announcements(owner, repo, category, count, token):
    """Return a list of dicts with keys: title, url, date, body_html."""

    query = """
    query($owner: String!, $repo: String!, $count: Int!, $category: String!) {
      repository(owner: $owner, name: $repo) {
        discussionCategories(first: 100) {
          nodes {
            name
            id
          }
        }
        discussions(
          first: $count,
          categoryId: null,
          orderBy: {field: CREATED_AT, direction: DESC}
        ) {
          nodes {
            title
            url
            createdAt
            bodyHTML
            category {
              name
            }
          }
        }
      }
    }
    """

    # First, we fetch discussions and filter by category name on the client
    # side, because the categoryId filter requires an ID we don't have yet.
    # An alternative is a two-step query; we keep it simple here by
    # over-fetching and filtering.
    #
    # We request more than `count` to have headroom after filtering.
    fetch_count = count * 10

    # Simpler single-step query: fetch recent discussions, filter locally
    query = """
    query($owner: String!, $repo: String!, $count: Int!) {
      repository(owner: $owner, name: $repo) {
        discussions(
          first: $count,
          orderBy: {field: CREATED_AT, direction: DESC}
        ) {
          nodes {
            title
            url
            createdAt
            bodyHTML
            category {
              name
            }
          }
        }
      }
    }
    """

    variables = {
        "owner": owner,
        "repo": repo,
        "count": fetch_count,
    }

    payload = json.dumps({"query": query, "variables": variables}).encode()

    req = urllib.request.Request(
        GITHUB_GRAPHQL_URL,
        data=payload,
        headers={
            "Authorization": f"bearer {token}",
            "Content-Type": "application/json",
        },
    )

    try:
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        body = exc.read().decode() if exc.fp else ""
        print(f"GitHub API error {exc.code}: {body}", file=sys.stderr)
        sys.exit(1)

    if "errors" in data:
        print(f"GraphQL errors: {json.dumps(data['errors'], indent=2)}",
              file=sys.stderr)
        sys.exit(1)

    discussions = data["data"]["repository"]["discussions"]["nodes"]

    # Filter to the requested category
    filtered = [
        d for d in discussions
        if d["category"]["name"].lower() == category.lower()
    ]

    # Take only the requested number
    filtered = filtered[:count]

    results = []
    for d in filtered:
        created = datetime.fromisoformat(
            d["createdAt"].replace("Z", "+00:00"))
        results.append({
            "title": d["title"],
            "url": d["url"],
            "date": created.strftime("%B %Y"),
            "date_iso": d["createdAt"],
            "body_html": d["bodyHTML"],
        })

    return results


def render_html(announcements):
    """Return an HTML snippet (no <html>/<body>) for Doxygen \\htmlinclude."""

    if not announcements:
        return "<!-- No announcements found -->\n"

    lines = []
    lines.append('<div class="announcements">')
    lines.append("<h2>Recent Announcements</h2>")
    lines.append("<p>Sourced from "
                 '<a href="https://github.com/pulseq/pulseq/discussions/'
                 'categories/announcements">GitHub Discussions</a>.</p>')

    for a in announcements:
        lines.append('<div class="announcement-entry" '
                     'style="margin-bottom:1.5em;">')
        lines.append(f'  <h3><a href="{a["url"]}">{a["title"]}</a></h3>')
        lines.append(f'  <p style="color:#666; font-size:0.9em;">'
                     f'{a["date"]}</p>')
        lines.append(f'  <div class="announcement-body">{a["body_html"]}'
                     f'</div>')
        lines.append("</div>")
        lines.append("<hr/>")

    # Remove trailing <hr/>
    if lines and lines[-1] == "<hr/>":
        lines.pop()

    lines.append("</div>")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Fetch GitHub Discussion announcements for the website.")
    parser.add_argument("-n", "--count", type=int, default=DEFAULT_COUNT,
                        help=f"Number of announcements (default: "
                             f"{DEFAULT_COUNT})")
    parser.add_argument("-o", "--output", default="announcements.html",
                        help="Output HTML file (default: announcements.html)")
    parser.add_argument("--owner", default=DEFAULT_REPO_OWNER,
                        help=f"Repository owner (default: {DEFAULT_REPO_OWNER})")
    parser.add_argument("--repo", default=DEFAULT_REPO_NAME,
                        help=f"Repository name (default: {DEFAULT_REPO_NAME})")
    parser.add_argument("--category", default=DEFAULT_CATEGORY,
                        help=f"Discussion category (default: {DEFAULT_CATEGORY})")
    args = parser.parse_args()

    token = os.environ.get("GITHUB_TOKEN", "")
    if not token:
        print("Error: GITHUB_TOKEN environment variable must be set.\n"
              "Create a token at https://github.com/settings/tokens with "
              "public_repo (classic) or Discussions read access "
              "(fine-grained).", file=sys.stderr)
        sys.exit(1)

    announcements = fetch_announcements(
        args.owner, args.repo, args.category, args.count, token)

    html = render_html(announcements)

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            args.output)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"Wrote {len(announcements)} announcement(s) to {out_path}")


if __name__ == "__main__":
    main()
