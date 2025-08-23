import React, { useEffect, useState } from "react";
import { Anchor, Group, Text } from "@mantine/core";
import { IconListSearch } from "@tabler/icons-react";
import classes from "./TableOfContents.module.css";

interface Heading {
  id: string;
  text: string;
  level: number;
}

export function TableOfContents() {
  const [headings, setHeadings] = useState<Heading[]>([]);
  const [activeId, setActiveId] = useState<string | null>(null);

  useEffect(() => {
    // Extract headings dynamically
    const selector = "h1, h2, h3";
    const elements = Array.from(
      document.querySelectorAll(selector),
    ) as HTMLElement[];

    const collectedHeadings = elements
      .filter((el) => el.id)
      .map((el) => ({
        id: el.id,
        text: el.innerText,
        level: parseInt(el.tagName[1]), // h1 -> 1, h2 -> 2, etc.
      }));

    setHeadings(collectedHeadings);
  }, []);

  useEffect(() => {
    if (headings.length === 0) return;

    const observer = new IntersectionObserver(
      (entries) => {
        const visibleHeadings = entries
          .filter((entry) => entry.isIntersecting)
          .sort((a, b) => a.boundingClientRect.top - b.boundingClientRect.top);

        if (visibleHeadings.length > 0) {
          setActiveId(visibleHeadings[0].target.id);
        }
      },
      {
        rootMargin: "0px 0px -70% 0px",
        threshold: 0.1,
      },
    );

    headings.forEach((heading) => {
      const el = document.getElementById(heading.id);
      if (el) observer.observe(el);
    });

    return () => observer.disconnect();
  }, [headings]);

  if (headings.length === 0) {
    return null; // No TOC if no headings
  }

  return (
    <>
      <Group mb="md">
        <IconListSearch size={18} />
        <Text style={{ fontVariant: "small-caps" }}>table of contents</Text>
      </Group>
      {headings.map((heading) => (
        <Anchor
          key={heading.id}
          href={`#${heading.id}`}
          className={`${classes.tocItem} ${activeId === heading.id ? classes.active : ""}`}
          style={{
            paddingLeft: `${heading.level * 20}px`,
            display: "block",
            marginBottom: "0.25rem",
          }}
          onClick={(e) => {
            e.preventDefault();
            document
              .getElementById(heading.id)
              ?.scrollIntoView({ behavior: "smooth" });
          }}
        >
          {heading.text}
        </Anchor>
      ))}
    </>
  );
}
